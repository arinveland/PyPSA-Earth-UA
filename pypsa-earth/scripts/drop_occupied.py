import pypsa
import geopandas as gpd
from pathlib import Path


def main(snakemake):
    network_in = Path(snakemake.input.network)
    gadm_path = Path(snakemake.input.gadm)
    network_out = Path(snakemake.output.network)

    # GADM IDs of occupied regions to remove 
    banned_gadm_ids = [
        "UA.15_1",  # Luhansk
        "UA.6_1",   # Donetsk
        "UA.26_1",  # Zaporizhzhia
        "UA.9_1",   # Kherson
        "UA.4_1",   # Crimea
        "UA.20_1",  # Sevastopol
    ]

    print(f"[drop_occupied] Reading network: {network_in}")
    n = pypsa.Network(network_in)

    print(f"[drop_occupied] Reading GADM shapes: {gadm_path}")
    gadm = gpd.read_file(gadm_path)

    if gadm.crs is not None and gadm.crs.to_epsg() != 4326:
        gadm = gadm.to_crs(epsg=4326)

    if "GADM_ID" not in gadm.columns:
        raise ValueError(
            f"[drop_occupied] 'GADM_ID' not in GADM file columns: {list(gadm.columns)}"
        )

    banned = gadm[gadm["GADM_ID"].isin(banned_gadm_ids)].copy()
    if banned.empty:
        raise ValueError(
            f"[drop_occupied] No polygons found for banned IDs: {banned_gadm_ids}"
        )

    print("[drop_occupied] Banned polygons:",
          ", ".join(sorted(banned["GADM_ID"].unique())))

    # build GeoDataFrame of all buses  
    buses_df = n.buses.copy()
    buses_gdf = gpd.GeoDataFrame(
        buses_df,
        geometry=gpd.points_from_xy(buses_df.x, buses_df.y),
        crs="EPSG:4326",
    )

    # spatial join: which buses lie inside banned oblasts
    buses_in_banned = gpd.sjoin(
        buses_gdf,
        banned[["geometry"]],
        how="inner",
        predicate="within",
    )

    buses_to_drop = buses_in_banned.index.unique().tolist()
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
            n.mremove(component[:-1].capitalize(), orphan_idx.tolist())

    # lines and links
    drop_orphan_components("lines", ["bus0", "bus1"])
    drop_orphan_components("links", ["bus0", "bus1", "bus2"])

    # one-bus components
    drop_orphan_components("generators", ["bus"])
    drop_orphan_components("loads", ["bus"])
    drop_orphan_components("stores", ["bus"])
    drop_orphan_components("storage_units", ["bus"])

    n.consistency_check()

        # ---------------------------------------------------------------------
    # YET TO BE THOUROUGHLY TESTED: Cleanup other components
    # ---------------------------------------------------------------------
    removed_h2 = {}

    # n.components keys are class names like "Bus", "Link", "Store" etc.
    for class_name, meta in n.components.items():
        list_name = meta.get("list_name")
        if not list_name or not hasattr(n, list_name):
            continue

        df = getattr(n, list_name)
        if df is None or df.empty:
            continue

        h2_ids = df.index[df.index.astype(str).str.contains("H2")]
        if len(h2_ids):
            n.mremove(class_name, h2_ids.tolist())
            removed_h2[class_name] = int(len(h2_ids))

    if removed_h2:
        print("[drop_occupied] Removed H2 components (by type):")
        for k in sorted(removed_h2):
            print(f"  - {k}: {removed_h2[k]}")
    else:
        print("[drop_occupied] No H2 components found to remove.")

    buses_str = n.buses.index.astype(str)

    parent_buses = buses_str[buses_str.str.fullmatch(r"\d+")]
    parent_buses_set = set(parent_buses)

    connected_by_grid = set()

    if hasattr(n, "lines") and not n.lines.empty:
        connected_by_grid |= set(n.lines.bus0.astype(str))
        connected_by_grid |= set(n.lines.bus1.astype(str))

    if hasattr(n, "transformers") and not n.transformers.empty:
        connected_by_grid |= set(n.transformers.bus0.astype(str))
        connected_by_grid |= set(n.transformers.bus1.astype(str))

    isolated_parents = sorted([b for b in parent_buses_set if b not in connected_by_grid])

    # remove the parent + its battery bus (if present)
    buses_to_drop_isolated = set(isolated_parents)
    buses_to_drop_isolated |= {f"{b} battery" for b in isolated_parents if f"{b} battery" in set(buses_str)}

    if buses_to_drop_isolated:
        print(f"[drop_occupied] Removing {len(isolated_parents)} isolated parent buses "
              f"(and {len(buses_to_drop_isolated) - len(isolated_parents)} battery counterparts).")
        n.mremove("Bus", sorted(buses_to_drop_isolated))
    else:
        print("[drop_occupied] No isolated parent buses found to remove.")

    drop_orphan_components("lines", ["bus0", "bus1"])
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
