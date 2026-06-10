from pathlib import Path
from typing import Dict

import geopandas as gpd
import pandas as pd
import pypsa
from shapely.geometry import LineString


EPS = 1e-9
CANDIDATE_SOLVERS = ["gurobi", "highs", "glpk"]


def _coerce_positive_float(name: str, value) -> float:
    out = float(value)
    if out <= 0:
        raise ValueError(f"{name} must be > 0, got {out}")
    return out


def _coerce_fraction(name: str, value) -> float:
    out = float(value)
    if out < -EPS or out > 1 + EPS:
        raise ValueError(f"{name} must be in [0, 1], got {out}")
    return max(0.0, min(1.0, out))


def _shift_to_modelling_year(ts: pd.Timestamp, modelling_year: int) -> pd.Timestamp:
    # Preserve month/day/time from weather reference year while switching to modelling year.
    delta_years = int(modelling_year) - int(ts.year)
    shifted = ts + pd.DateOffset(years=delta_years)
    return pd.Timestamp(shifted)


def _load_cfg() -> Dict[str, object]:
    cfg = snakemake.config.get("scenario_analysis", {})
    if not isinstance(cfg, dict):
        raise ValueError("scenario_analysis config section must be a mapping")

    run_dir_name = str(cfg.get("rdir_name", cfg.get("run_name", ""))).strip()
    modelling_year = int(cfg.get("modelling_year", cfg.get("Modelling_year")))
    snapshots_start = pd.Timestamp(cfg.get("snapshots_start"))
    snapshots_end = pd.Timestamp(cfg.get("snapshots_end"))
    if snapshots_end <= snapshots_start:
        raise ValueError("snapshots_end must be strictly after snapshots_start")

    attacks_file = str(cfg.get("attacks_file", "")).strip()
    if not attacks_file:
        raise ValueError("scenario_analysis.attacks_file must be non-empty")

    gen_cfg = cfg.get("generation", {})
    trans_cfg = cfg.get("transmission", {})

    gen_radius_km = _coerce_positive_float("generation.gen_radius", gen_cfg.get("gen_radius"))
    gen_damage = _coerce_fraction(
        "generation.gen_damage_percentage", gen_cfg.get("gen_damage_percentage")
    )
    gen_days_to_repair = int(gen_cfg.get("days_to_repair"))
    if gen_days_to_repair < 0:
        raise ValueError("generation.days_to_repair must be >= 0")

    trans_radius_km = _coerce_positive_float("transmission.trans_radius", trans_cfg.get("trans_radius"))
    trans_damage = _coerce_fraction(
        "transmission.trans_damage_percentage", trans_cfg.get("trans_damage_percentage")
    )
    trans_days_to_repair = int(trans_cfg.get("days_to_repair"))
    if trans_days_to_repair < 0:
        raise ValueError("transmission.days_to_repair must be >= 0")

    line_capacity_factor = float(snakemake.config.get("line_capacity_factor", 0.7))

    return {
        "run_dir_name": run_dir_name,
        "modelling_year": modelling_year,
        "snapshots_start": snapshots_start,
        "snapshots_end": snapshots_end,
        "attacks_file": attacks_file,
        "gen_radius_km": gen_radius_km,
        "gen_damage": gen_damage,
        "gen_days_to_repair": gen_days_to_repair,
        "trans_radius_km": trans_radius_km,
        "trans_damage": trans_damage,
        "trans_days_to_repair": trans_days_to_repair,
        "line_capacity_factor": line_capacity_factor,
    }


def _apply_snapshot_window(n: pypsa.Network, start: pd.Timestamp, end: pd.Timestamp) -> None:
    snaps = pd.DatetimeIndex(n.snapshots)
    selected = snaps[(snaps >= start) & (snaps < end)]
    if selected.empty:
        raise ValueError("No snapshots remain after applying [snapshots_start, snapshots_end)")
    n.set_snapshots(selected)


def _read_attacks(attacks_path: Path) -> pd.DataFrame:
    if not attacks_path.exists():
        raise FileNotFoundError(f"Attacks file not found: {attacks_path}")

    df = pd.read_csv(attacks_path)
    required = {"EVENT_DATE", "LATITUDE", "LONGITUDE"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Attacks file missing columns: {sorted(missing)}")

    df = df.copy()
    df["EVENT_DATE"] = pd.to_datetime(df["EVENT_DATE"], errors="coerce")
    df = df.dropna(subset=["EVENT_DATE", "LATITUDE", "LONGITUDE"])
    if df.empty:
        raise ValueError("No valid attack rows after parsing EVENT_DATE/LATITUDE/LONGITUDE")

    return df


def _filter_attacks_window(
    attacks: pd.DataFrame,
    modelling_start: pd.Timestamp,
    days_to_repair: int,
) -> pd.DataFrame:
    window_start = modelling_start - pd.Timedelta(days=int(days_to_repair))
    window_end = modelling_start + pd.Timedelta(days=1)
    return attacks[(attacks["EVENT_DATE"] >= window_start) & (attacks["EVENT_DATE"] < window_end)].copy()


def _attacks_gdf(attacks: pd.DataFrame) -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        attacks.copy(),
        geometry=gpd.points_from_xy(attacks["LONGITUDE"], attacks["LATITUDE"]),
        crs="EPSG:4326",
    )


def _generation_hit_counts(
    n: pypsa.Network,
    attacks_gdf: gpd.GeoDataFrame,
    radius_km: float,
) -> pd.Series:
    gen_buses = set(n.generators["bus"].unique()) if not n.generators.empty else set()

    storage_buses: set = set()
    if not n.storage_units.empty and "carrier" in n.storage_units.columns:
        carriers = n.storage_units["carrier"].astype(str).str.strip().str.lower()
        storage_buses = set(n.storage_units.loc[carriers.isin({"hydro", "phs"}), "bus"].unique())

    all_buses = gen_buses | storage_buses
    if not all_buses:
        return pd.Series(dtype=int)

    bus_points = n.buses.loc[n.buses.index.isin(all_buses), ["x", "y"]].copy()
    if bus_points.empty:
        return pd.Series(dtype=int)

    bus_gdf = gpd.GeoDataFrame(
        bus_points,
        geometry=gpd.points_from_xy(bus_points["x"], bus_points["y"]),
        crs="EPSG:4326",
    )

    attacks_m = attacks_gdf.to_crs(epsg=3857)
    buses_m = bus_gdf.to_crs(epsg=3857)

    buffered = attacks_m[["geometry"]].copy()
    buffered["attack_id"] = attacks_m.index.astype(str)
    buffered["geometry"] = buffered.geometry.buffer(radius_km * 1000.0)

    hits = gpd.sjoin(
        buses_m,
        buffered[["attack_id", "geometry"]],
        how="inner",
        predicate="within",
    )
    if hits.empty:
        return pd.Series(dtype=int)

    return hits.groupby(hits.index)["attack_id"].nunique().astype(int)


def _line_like_geometries(n: pypsa.Network, component: str) -> gpd.GeoDataFrame:
    df = getattr(n, component)
    if df.empty:
        return gpd.GeoDataFrame(columns=["geometry"], geometry="geometry", crs="EPSG:4326")

    rows = []
    for asset_id, row in df.iterrows():
        bus0 = row.get("bus0")
        bus1 = row.get("bus1")
        if bus0 not in n.buses.index or bus1 not in n.buses.index:
            continue

        x0 = float(n.buses.at[bus0, "x"])
        y0 = float(n.buses.at[bus0, "y"])
        x1 = float(n.buses.at[bus1, "x"])
        y1 = float(n.buses.at[bus1, "y"])
        geom = LineString([(x0, y0), (x1, y1)])
        rows.append((asset_id, geom))

    if not rows:
        return gpd.GeoDataFrame(columns=["geometry"], geometry="geometry", crs="EPSG:4326")

    out = gpd.GeoDataFrame(rows, columns=["asset_id", "geometry"], geometry="geometry", crs="EPSG:4326")
    out = out.set_index("asset_id")
    return out


def _transmission_hit_counts(
    n: pypsa.Network,
    attacks_gdf: gpd.GeoDataFrame,
    radius_km: float,
    component: str,
) -> pd.Series:
    assets = _line_like_geometries(n, component)
    if assets.empty:
        return pd.Series(dtype=int)

    attacks_m = attacks_gdf.to_crs(epsg=3857)
    assets_m = assets.to_crs(epsg=3857)

    buffered = attacks_m[["geometry"]].copy()
    buffered["attack_id"] = attacks_m.index.astype(str)
    buffered["geometry"] = buffered.geometry.buffer(radius_km * 1000.0)

    hits = gpd.sjoin(
        assets_m,
        buffered[["attack_id", "geometry"]],
        how="inner",
        predicate="intersects",
    )
    if hits.empty:
        return pd.Series(dtype=int)

    return hits.groupby(hits.index)["attack_id"].nunique().astype(int)


def _apply_damage_with_counts(series: pd.Series, counts: pd.Series, damage_fraction: float) -> pd.Series:
    if series.empty or counts.empty or damage_fraction <= EPS:
        return series

    out = series.copy().astype(float)
    common = out.index.intersection(counts.index)
    if common.empty:
        return out

    factors = (1.0 - float(damage_fraction)) ** counts.loc[common].astype(int)
    out.loc[common] = out.loc[common] * factors
    out.loc[out < 0.0] = 0.0
    return out


def _set_line_capacity_factor(n: pypsa.Network, factor: float) -> None:
    if n.lines.empty or "s_max_pu" not in n.lines.columns:
        return
    n.lines.loc[:, "s_max_pu"] = factor
    print(f"[scenario_analysis] Set s_max_pu = {factor} on {len(n.lines)} lines")


def _freeze_expansion(n: pypsa.Network) -> None:
    for comp_name, flag_col in [
        ("generators", "p_nom_extendable"),
        ("storage_units", "p_nom_extendable"),
        ("links", "p_nom_extendable"),
        ("stores", "e_nom_extendable"),
        ("lines", "s_nom_extendable"),
        ("transformers", "s_nom_extendable"),
    ]:
        comp_df = getattr(n, comp_name, None)
        if comp_df is not None and not comp_df.empty and flag_col in comp_df.columns:
            comp_df.loc[:, flag_col] = False


def _set_storage_cyclic(n: pypsa.Network) -> None:
    if not n.storage_units.empty and "cyclic_state_of_charge" in n.storage_units.columns:
        n.storage_units.loc[:, "cyclic_state_of_charge"] = True


def _solve_network(n: pypsa.Network) -> str:
    last_error = None
    for solver in CANDIDATE_SOLVERS:
        try:
            if hasattr(n, "optimize"):
                # `extra_functionality` may be attached to the Network before calling this
                # function by setting `n._extra_functionality` on the network object.
                ef = getattr(n, "_extra_functionality", None)
                result = n.optimize(solver_name=solver, extra_functionality=ef)
                if isinstance(result, tuple) and len(result) >= 2:
                    status = str(result[0]).lower()
                    condition = str(result[1]).lower()
                    if "ok" in status or "optimal" in condition:
                        return solver
                else:
                    return solver
            elif hasattr(n, "lopf"):
                ef = getattr(n, "_extra_functionality", None)
                ok = n.lopf(pyomo=False, solver_name=solver, extra_functionality=ef)
                if ok is not False:
                    return solver
            else:
                raise RuntimeError("No supported PyPSA optimization method found")
        except Exception as err:  # noqa: BLE001
            last_error = err
            continue

    raise RuntimeError(f"Optimization failed for all candidate solvers: {last_error}")


def main() -> None:
    cfg = _load_cfg()

    n = pypsa.Network(snakemake.input["network"])
    _apply_snapshot_window(n, cfg["snapshots_start"], cfg["snapshots_end"])

    modelling_start = _shift_to_modelling_year(cfg["snapshots_start"], cfg["modelling_year"])

    attacks_path = Path("Files") / cfg["attacks_file"]
    attacks = _read_attacks(attacks_path)

    gen_attacks = _filter_attacks_window(
        attacks,
        modelling_start,
        cfg["gen_days_to_repair"],
    )
    trans_attacks = _filter_attacks_window(
        attacks,
        modelling_start,
        cfg["trans_days_to_repair"],
    )

    gen_hits = pd.Series(dtype=int)
    if not gen_attacks.empty:
        gen_hits = _generation_hit_counts(
            n,
            _attacks_gdf(gen_attacks),
            cfg["gen_radius_km"],
        )

    trans_hits_lines = pd.Series(dtype=int)
    trans_hits_transformers = pd.Series(dtype=int)
    if not trans_attacks.empty:
        trans_gdf = _attacks_gdf(trans_attacks)
        trans_hits_lines = _transmission_hit_counts(
            n,
            trans_gdf,
            cfg["trans_radius_km"],
            component="lines",
        )
        trans_hits_transformers = _transmission_hit_counts(
            n,
            trans_gdf,
            cfg["trans_radius_km"],
            component="transformers",
        )

    if not n.generators.empty:
        bus_hit_counts = gen_hits if not gen_hits.empty else pd.Series(dtype=int)
        generator_counts = n.generators["bus"].map(bus_hit_counts).fillna(0).astype(int)
        non_nuclear_mask = pd.Series(True, index=n.generators.index)
        if "carrier" in n.generators.columns:
            # Nuclear and cross-border units are assumed operational and excluded from damage.
            carriers = n.generators["carrier"].astype(str).str.strip().str.lower()
            non_nuclear_mask = ~carriers.isin({"nuclear", "cross-border"})

        damage_targets = generator_counts[(generator_counts > 0) & non_nuclear_mask]
        n.generators.loc[:, "p_nom"] = _apply_damage_with_counts(
            n.generators["p_nom"],
            damage_targets,
            cfg["gen_damage"],
        )

    if not n.storage_units.empty:
        bus_hit_counts = gen_hits if not gen_hits.empty else pd.Series(dtype=int)
        storage_counts = n.storage_units["bus"].map(bus_hit_counts).fillna(0).astype(int)
        hydro_mask = pd.Series(False, index=n.storage_units.index)
        if "carrier" in n.storage_units.columns:
            carriers = n.storage_units["carrier"].astype(str).str.strip().str.lower()
            hydro_mask = carriers.isin({"hydro", "phs"})
        damage_targets_storage = storage_counts[(storage_counts > 0) & hydro_mask]
        n.storage_units.loc[:, "p_nom"] = _apply_damage_with_counts(
            n.storage_units["p_nom"],
            damage_targets_storage,
            cfg["gen_damage"],
        )

    for comp in ["lines", "transformers"]:
        hits = trans_hits_lines if comp == "lines" else trans_hits_transformers
        df = getattr(n, comp)
        if df.empty:
            continue
        if "s_nom" in df.columns:
            df.loc[:, "s_nom"] = _apply_damage_with_counts(df["s_nom"], hits, cfg["trans_damage"])
        if "s_nom_min" in df.columns:
            df.loc[:, "s_nom_min"] = _apply_damage_with_counts(
                df["s_nom_min"], hits, cfg["trans_damage"]
            )
        if "s_nom_max" in df.columns:
            df.loc[:, "s_nom_max"] = _apply_damage_with_counts(
                df["s_nom_max"], hits, cfg["trans_damage"]
            )

    _freeze_expansion(n)
    _set_storage_cyclic(n)
    _set_line_capacity_factor(n, cfg["line_capacity_factor"])
    added_ls = n.optimize.add_load_shedding(sign=1, marginal_cost=100_000.0, p_nom=1e6)
    if len(added_ls):
        n.generators.loc[added_ls, "carrier"] = "load shedding"
        if "load shedding" not in n.carriers.index:
            n.add("Carrier", "load shedding", color="#c20808")
        print(f"[scenario_analysis] Added {len(added_ls)} load-shedding generators")
    # Attach an extra_functionality hook to enforce cross-border import limits
    def _add_cross_border_limit(network: pypsa.Network, snapshots) -> None:
        cb_cfg = snakemake.config.get("cross_border", {})
        if not cb_cfg.get("enable", False):
            return

        carrier = cb_cfg.get("carrier", "cross-border")
        ntc_gw = float(cb_cfg.get("ntc_gw", 0.0))
        ntc_mw = ntc_gw * 1000.0

        # Match carrier strings case-insensitively
        gens_mask = (
            network.generators["carrier"].astype(str).str.strip().str.lower()
            == str(carrier).strip().lower()
        )
        gens = network.generators.index[gens_mask]

        if len(gens) == 0:
            print(
                f"[extra_functionality] cross_border.enable is True but no generators "
                f"with carrier {carrier!r} were found; skipping cross-border limit."
            )
            return

        try:
            m = network.model
            p_gen = m.variables["Generator-p"].sel(snapshot=snapshots, Generator=gens)
            expr = p_gen.sum("Generator")

            # Bidirectional capacity limit per snapshot
            m.add_constraints(expr <= ntc_mw, name="max_cross_border_plus")
            m.add_constraints(expr >= -ntc_mw, name="max_cross_border_minus")

            print(
                f"[extra_functionality] Added per-snapshot cross-border limit ±{ntc_mw:.1f} MW "
                f"for carrier {carrier!r}."
            )
        except Exception as exc:  # noqa: BLE001
            print(f"[extra_functionality] Failed to add cross-border constraint: {exc}")

    # store hook on network so _solve_network can pick it up
    n._extra_functionality = _add_cross_border_limit
    _solve_network(n)

    solved_network_path = Path(snakemake.output["network"])
    solved_network_path.parent.mkdir(parents=True, exist_ok=True)

    expected_dir = Path("Files") / "results"
    if cfg["run_dir_name"]:
        expected_dir = expected_dir / cfg["run_dir_name"]

    expected_network_path = expected_dir / "scenario_analysis.nc"
    if solved_network_path.resolve() != expected_network_path.resolve():
        raise RuntimeError(
            f"Snakefile output path {solved_network_path} does not match expected path {expected_network_path}. "
            "Update the Snakefile so ua_scenario_analysis writes to Files/results/<rdir_name>/scenario_analysis.nc."
        )

    n.export_to_netcdf(solved_network_path)
    print(f"[scenario_analysis] Exported solved network to {solved_network_path}")


if __name__ == "__main__":
    main()
