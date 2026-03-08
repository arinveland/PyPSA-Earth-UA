import pypsa
import geopandas as gpd
from pathlib import Path


def main(snakemake):
    network_in = Path(snakemake.input.network)
    gadm_path = Path(snakemake.input.gadm)
    network_out = Path(snakemake.output.network)

    cfg = snakemake.config.get("drop_occupied", {})
    if not cfg.get("enable", False):
        n.export_to_netcdf(snakemake.output["network"])
        raise SystemExit()

    # GADM IDs of occupied regions to remove 
    occupied_ids = cfg.get("occupied_ids", [])

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

    # write out pruned network 
    network_out.parent.mkdir(parents=True, exist_ok=True)
    print(f"[drop_occupied] Writing pruned network to: {network_out}")
    n.export_to_netcdf(network_out)


if __name__ == "__main__":
    main(snakemake)
