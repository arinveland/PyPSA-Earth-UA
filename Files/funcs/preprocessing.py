from pathlib import Path
import pandas as pd

def preprocess(network, original_path):
    n = network.copy()

    #Remove H2 components
    removed_counts = {}

    for class_name, meta in n.components.items():
        list_name = meta.get("list_name")
        if not list_name or not hasattr(n, list_name):
            continue
    
        df = getattr(n, list_name)
        if df is None or df.empty:
            continue
        
        h2_ids = df.index[df.index.astype(str).str.contains("H2")]
        if len(h2_ids) == 0:
            continue
    
        for name in list(map(str, h2_ids)):
            n.remove(class_name, name)
    
    removed_counts[class_name] = int(len(h2_ids))

    #Remove isolated nodes
    parent_buses = n.buses.index.astype(str)
    parent_buses = parent_buses[parent_buses.str.fullmatch(r"\d+")]

    connected_by_lines = set()

    if hasattr(n, "lines") and not n.lines.empty:
        connected_by_lines |= set(n.lines["bus0"].astype(str))
        connected_by_lines |= set(n.lines["bus1"].astype(str))

    if hasattr(n, "transformers") and not n.transformers.empty:
        connected_by_lines |= set(n.transformers["bus0"].astype(str))
        connected_by_lines |= set(n.transformers["bus1"].astype(str))

    isolated_parents = [b for b in parent_buses if b not in connected_by_lines]

    to_remove_buses = set(isolated_parents)
    to_remove_buses |= {f"{b} battery" for b in isolated_parents if f"{b} battery" in n.buses.index.astype(str)}

    removed_attached = remove_components(n, to_remove_buses)

    removed_bus_count = 0
    for b in list(to_remove_buses):
        if b in n.buses.index.astype(str).tolist():
            n.remove("Bus", b)
            removed_bus_count += 1

    new_path = original_path.removesuffix(".nc") + "_preprocessed.nc"
    n.export_to_netcdf(new_path)

    print("Removed components (by type):")
    for k in sorted(removed_counts):
        print(f"  {k}: {removed_counts[k]}")
    print(f"Isolated parent buses found (no Lines/Transformers): {len(isolated_parents)}")
    print(f"Total buses removed (parents + battery counterparts): {removed_bus_count}")
    print("Removed attached components (by type):")
    for k in sorted(removed_attached):
        print(f"  {k}: {removed_attached[k]}")
    print(f"\nSaved pruned network to:\n{new_path}")


#Helper function to remove any components attached to a bus
def remove_components(nw, buses_to_remove):
    buses_to_remove = set(map(str, buses_to_remove))
    removed = {}

    # iterate component classes
    for class_name, meta in nw.components.items():
        list_name = meta.get("list_name")
        if not list_name or not hasattr(nw, list_name):
            continue

        df = getattr(nw, list_name)
        if df is None or df.empty:
            continue

        # remove components that reference these buses in bus-like columns
        bus_cols = [c for c in ["bus", "bus0", "bus1", "bus2", "bus3"] if c in df.columns]
        if not bus_cols:
            continue

        mask = pd.Series(False, index=df.index)
        for c in bus_cols:
            mask |= df[c].astype(str).isin(buses_to_remove)

        ids = df.index[mask]
        if len(ids) == 0:
            continue

        # Network.remove in your version removes one-by-one (expects single name)
        for name in map(str, ids):
            nw.remove(class_name, name)

        removed[class_name] = int(len(ids))

    return removed
