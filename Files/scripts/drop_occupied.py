import pypsa
import geopandas as gpd
from pathlib import Path
import math


def _line_transformer_components(n):
    # Build connected components using only lines and transformers.
    adjacency = {bus: set() for bus in n.buses.index}

    for comp_name in ["lines", "transformers"]:
        df = getattr(n, comp_name)
        if df.empty:
            continue
        for _, row in df.iterrows():
            bus0 = row.get("bus0")
            bus1 = row.get("bus1")
            if bus0 in adjacency and bus1 in adjacency:
                adjacency[bus0].add(bus1)
                adjacency[bus1].add(bus0)

    seen = set()
    components = []
    for root in adjacency:
        if root in seen:
            continue
        stack = [root]
        comp = set()
        seen.add(root)
        while stack:
            node = stack.pop()
            comp.add(node)
            for nbr in adjacency[node]:
                if nbr not in seen:
                    seen.add(nbr)
                    stack.append(nbr)
        components.append(comp)

    return components


def _nearest_main_bus_map(n, island_buses, main_buses):
    # Map each island bus to the closest bus in the main connected component.
    main_xy = n.buses.loc[list(main_buses), ["x", "y"]]
    island_xy = n.buses.loc[list(island_buses), ["x", "y"]]
    if main_xy.empty or island_xy.empty:
        return {}

    mapping = {}
    for island_bus, row in island_xy.iterrows():
        ix = float(row["x"])
        iy = float(row["y"])

        # Euclidean distance in lon/lat space is sufficient for nearest-bus selection here.
        distances = (main_xy["x"] - ix) ** 2 + (main_xy["y"] - iy) ** 2
        nearest = distances.idxmin()
        if isinstance(nearest, float) and math.isnan(nearest):
            continue
        mapping[island_bus] = nearest

    return mapping


def _reassign_components_from_islands(n, bus_mapping):
    if not bus_mapping:
        return

    # Any component bus reference pointing at an island bus is reassigned.
    component_bus_cols = {
        "lines": ["bus0", "bus1"],
        "transformers": ["bus0", "bus1"],
        "links": ["bus0", "bus1", "bus2"],
        "generators": ["bus"],
        "loads": ["bus"],
        "stores": ["bus"],
        "storage_units": ["bus"],
    }

    for component, bus_cols in component_bus_cols.items():
        df = getattr(n, component)
        if df.empty:
            continue
        for col in bus_cols:
            if col not in df.columns:
                continue
            mask = df[col].isin(bus_mapping)
            if mask.any():
                df.loc[mask, col] = df.loc[mask, col].map(bus_mapping)


def _drop_degenerate_branches(n):
    # After reassignment, branches may collapse to the same endpoint bus.
    for component, bus0_col, bus1_col, pypsa_name in [
        ("lines", "bus0", "bus1", "Line"),
        ("transformers", "bus0", "bus1", "Transformer"),
    ]:
        df = getattr(n, component)
        if df.empty or bus0_col not in df.columns or bus1_col not in df.columns:
            continue
        degenerate = df.index[df[bus0_col] == df[bus1_col]]
        if len(degenerate):
            print(f"[drop_occupied] Dropping {len(degenerate)} degenerate {component} after reassignment")
            n.mremove(pypsa_name, degenerate.tolist())


def main(snakemake):
    network_in = Path(snakemake.input.network)
    gadm_path = Path(snakemake.input.gadm)
    network_out = Path(snakemake.output.network)
    n = pypsa.Network(network_in)

    cfg = snakemake.config.get("drop_occupied", {})
    if not cfg.get("enable", False):
        n.export_to_netcdf(network_out)
        raise SystemExit()

    # GADM IDs of occupied regions to remove 
    occupied_ids = cfg.get("occupied_ids", [])

    print(f"[drop_occupied] Reading network: {network_in}")

    print(f"[drop_occupied] Reading GADM shapes: {gadm_path}")
    gadm = gpd.read_file(gadm_path)

    if gadm.crs is not None and gadm.crs.to_epsg() != 4326:
        gadm = gadm.to_crs(epsg=4326)

    if "GADM_ID" not in gadm.columns:
        raise ValueError(
            f"[drop_occupied] 'GADM_ID' not in GADM file columns: {list(gadm.columns)}"
        )

    occupied = gadm[gadm["GADM_ID"].isin(occupied_ids)].copy()
    if occupied.empty:
        raise ValueError(
            f"[drop_occupied] No polygons found for occupied IDs: {occupied_ids}"
        )

    print("[drop_occupied] Occupied polygons:",
          ", ".join(sorted(occupied["GADM_ID"].unique())))

    # build GeoDataFrame of all buses  
    buses_df = n.buses.copy()
    buses_gdf = gpd.GeoDataFrame(
        buses_df,
        geometry=gpd.points_from_xy(buses_df.x, buses_df.y),
        crs="EPSG:4326",
    )

    # spatial join: which buses lie inside occupied oblasts
    buses_in_occupied = gpd.sjoin(
        buses_gdf,
        occupied[["geometry"]],
        how="inner",
        predicate="within",
    )

    buses_to_drop = buses_in_occupied.index.unique().tolist()
    print(f"[drop_occupied] Found {len(buses_to_drop)} buses to drop.")

    # remove relevant buses (cascade-deletes attached components) 
    if buses_to_drop:
        n.mremove("Bus", buses_to_drop)
        print(
            f"[drop_occupied] After n.mremove('Bus', ...): "
            f"{len(n.buses)} buses, {len(n.lines)} lines, "
            f"{len(n.links)} links, {len(n.generators)} generators"
        )
    else:
        print("[drop_occupied] No buses removed.")

    # cleanup of orphan components

    # helper to drop components whose bus columns point to non-existent buses
    def drop_orphan_components(component, bus_cols):
        df = getattr(n, component)
        if df.empty:
            return
        orphan_mask = False
        for col in bus_cols:
            if col in df.columns:
                orphan_mask |= ~df[col].isin(n.buses.index)
        orphan_idx = df.index[orphan_mask]
        if len(orphan_idx):
            print(f"[drop_occupied] Dropping {len(orphan_idx)} {component} with invalid buses")
            component_map = {
                "lines": "Line",
                "transformers": "Transformer",
                "links": "Link",
                "generators": "Generator",
                "loads": "Load",
                "stores": "Store",
                "storage_units": "StorageUnit",
            }
            n.mremove(component_map[component], orphan_idx.tolist())

    # lines and links
    drop_orphan_components("lines", ["bus0", "bus1"])
    drop_orphan_components("transformers", ["bus0", "bus1"])
    drop_orphan_components("links", ["bus0", "bus1", "bus2"])

    # one-bus components
    drop_orphan_components("generators", ["bus"])
    drop_orphan_components("loads", ["bus"])
    drop_orphan_components("stores", ["bus"])
    drop_orphan_components("storage_units", ["bus"])

    # Rewire buses in disconnected line/transformer islands to nearest bus in main network,
    # then remove the island buses.
    components = _line_transformer_components(n)
    if len(components) > 1:
        main_component = max(components, key=len)
        island_buses = set().union(*[c for c in components if c is not main_component])
        print(
            f"[drop_occupied] Found {len(components) - 1} island(s) with "
            f"{len(island_buses)} bus(es); rewiring to main network."
        )

        bus_mapping = _nearest_main_bus_map(n, island_buses, main_component)
        _reassign_components_from_islands(n, bus_mapping)

        buses_to_remove = sorted(set(bus_mapping.keys()).intersection(n.buses.index))
        if buses_to_remove:
            n.mremove("Bus", buses_to_remove)
            print(f"[drop_occupied] Removed {len(buses_to_remove)} island bus(es) after reassignment")

        _drop_degenerate_branches(n)

        # Final orphan cleanup in case reassignment/removal invalidated any references.
        drop_orphan_components("lines", ["bus0", "bus1"])
        drop_orphan_components("transformers", ["bus0", "bus1"])
        drop_orphan_components("links", ["bus0", "bus1", "bus2"])
        drop_orphan_components("generators", ["bus"])
        drop_orphan_components("loads", ["bus"])
        drop_orphan_components("stores", ["bus"])
        drop_orphan_components("storage_units", ["bus"])

    n.consistency_check()

    # write out pruned network 
    network_out.parent.mkdir(parents=True, exist_ok=True)
    print(f"[drop_occupied] Writing pruned network to: {network_out}")
    n.export_to_netcdf(network_out)


if __name__ == "__main__":
    main(snakemake)
