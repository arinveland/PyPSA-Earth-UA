from pathlib import Path

import geopandas as gpd
import pypsa


EPS = 1e-9


def _coerce_mapping(name, value):
	if isinstance(value, dict):
		return value
	if isinstance(value, list):
		merged = {}
		for item in value:
			if not isinstance(item, dict):
				raise ValueError(
					f"{name} must be a mapping or list of single-item mappings, got item {item!r}"
				)
			for key, val in item.items():
				merged[str(key)] = val
		return merged
	raise ValueError(f"{name} must be a mapping, got {type(value).__name__}")


def _bus_region_mapping(network, gadm_path, configured_regions):
	gadm = gpd.read_file(gadm_path)
	if "GADM_ID" not in gadm.columns:
		raise ValueError(f"'GADM_ID' column not found in {gadm_path}")

	gadm = gadm[gadm["GADM_ID"].isin(configured_regions)].copy()
	if gadm.empty:
		raise ValueError("No configured load_adjustment regions found in GADM shapes")

	if gadm.crs is not None and gadm.crs.to_epsg() != 4326:
		gadm = gadm.to_crs(epsg=4326)

	buses_df = network.buses.copy()
	if "x" not in buses_df.columns or "y" not in buses_df.columns:
		raise ValueError("Network buses must include x/y coordinates")

	buses_gdf = gpd.GeoDataFrame(
		buses_df,
		geometry=gpd.points_from_xy(buses_df.x, buses_df.y),
		crs="EPSG:4326",
	)

	joined = gpd.sjoin(
		buses_gdf,
		gadm[["GADM_ID", "geometry"]],
		how="left",
		predicate="intersects",
	)

	joined = joined.copy()
	joined["bus_id"] = joined.index
	joined = joined.sort_values(["bus_id", "GADM_ID"], na_position="last")
	joined = joined.drop_duplicates(subset=["bus_id"], keep="first")

	return dict(zip(joined["bus_id"], joined["GADM_ID"]))


def _coerce_numeric_mapping(cfg, key_name):
	raw = _coerce_mapping(f"load_adjustment.{key_name}", cfg.get(key_name, {}))
	out = {}
	for region, value in raw.items():
		try:
			out[str(region)] = float(value)
		except (TypeError, ValueError) as exc:
			raise ValueError(
				f"load_adjustment.{key_name}.{region} must be numeric, got {value!r}"
			) from exc
	return out


def main(snakemake):
	network_in = Path(snakemake.input["network"])
	gadm_path = Path(snakemake.input["gadm"])
	network_out = Path(snakemake.output["network"])

	n = pypsa.Network(network_in)

	cfg = snakemake.config.get("load_adjustment", {})
	if not cfg.get("enable", False):
		n.export_to_netcdf(network_out)
		raise SystemExit()

	if n.loads.empty:
		n.export_to_netcdf(network_out)
		raise SystemExit()

	base_industrial = _coerce_numeric_mapping(cfg, "base_industrial_loads")
	industrial_coeff = _coerce_numeric_mapping(cfg, "industrial_adjustment_coefficients")
	household_coeff = _coerce_numeric_mapping(cfg, "population_adjustment_coefficients")

	configured_regions = set(base_industrial) | set(industrial_coeff) | set(household_coeff)
	if not configured_regions:
		raise ValueError("load_adjustment config is enabled, but no region mappings were provided")

	missing_base = configured_regions - set(base_industrial)
	missing_ind = configured_regions - set(industrial_coeff)
	missing_house = configured_regions - set(household_coeff)

	if missing_base:
		raise ValueError(
			"Missing load_adjustment.base_industrial_loads entries for regions: "
			f"{sorted(missing_base)}"
		)
	if missing_ind:
		raise ValueError(
			"Missing load_adjustment.industrial_adjustment_coefficients entries for regions: "
			f"{sorted(missing_ind)}"
		)
	if missing_house:
		raise ValueError(
			"Missing load_adjustment.population_adjustment_coefficients entries for regions: "
			f"{sorted(missing_house)}"
		)

	for region in sorted(configured_regions):
		if base_industrial[region] < -EPS:
			raise ValueError(
				"load_adjustment.base_industrial_loads values must be >= 0; "
				f"got {base_industrial[region]} for region {region}"
			)
		if industrial_coeff[region] < -EPS:
			raise ValueError(
				"load_adjustment.industrial_adjustment_coefficients values must be >= 0; "
				f"got {industrial_coeff[region]} for region {region}"
			)
		if household_coeff[region] < -EPS:
			raise ValueError(
				"load_adjustment.population_adjustment_coefficients values must be >= 0; "
				f"got {household_coeff[region]} for region {region}"
			)

	bus_to_region = _bus_region_mapping(
		network=n,
		gadm_path=gadm_path,
		configured_regions=configured_regions,
	)

	loads = n.loads.copy()
	load_region = loads["bus"].map(bus_to_region)
	load_region = load_region.where(load_region.notna(), None)

	region_to_loads = {region: [] for region in configured_regions}
	for load_name, region in load_region.items():
		if region in region_to_loads:
			region_to_loads[region].append(load_name)

	p_set = n.loads_t.p_set if hasattr(n.loads_t, "p_set") else None
	has_ts = p_set is not None and not p_set.empty

	region_scale = {}
	for region in sorted(configured_regions):
		load_names = region_to_loads[region]
		if not load_names:
			raise ValueError(
				f"No loads found in configured region {region}. "
				"Update load_adjustment config or region mapping."
			)

		if has_ts:
			existing_total = float(p_set[load_names].sum().sum())
		else:
			existing_total = float(loads.loc[load_names, "p_set"].sum())

		industrial_base = base_industrial[region] * 1e6
		if industrial_base > existing_total + 1e-6:
			industrial_base = 0.70 * max(0.0, existing_total)

		household_base = max(0.0, existing_total - industrial_base)
		adjusted_total = (
			industrial_base * industrial_coeff[region]
			+ household_base * household_coeff[region]
		)

		if existing_total <= EPS:
			if adjusted_total > 1e-6:
				raise ValueError(
					f"Region {region} has zero current load but positive adjusted load target"
				)
			region_scale[region] = 1.0
		else:
			region_scale[region] = adjusted_total / existing_total

	for region, factor in region_scale.items():
		load_names = region_to_loads[region]
		loads.loc[load_names, "p_set"] = loads.loc[load_names, "p_set"] * factor
		if has_ts:
			p_set.loc[:, load_names] = p_set.loc[:, load_names] * factor

	n.loads["p_set"] = loads["p_set"]
	if has_ts:
		n.loads_t.p_set = p_set

	n.consistency_check()
	network_out.parent.mkdir(parents=True, exist_ok=True)
	n.export_to_netcdf(network_out)


if __name__ == "__main__":
	main(snakemake)  # type: ignore[name-defined]
