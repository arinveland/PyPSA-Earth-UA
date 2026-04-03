from collections import defaultdict, deque

import geopandas as gpd
import pypsa


EPS = 1e-9


def _normalize_token(value):
	return "".join(ch for ch in str(value).lower() if ch.isalnum())


def _nearly_equal(a, b, tol=1e-6):
	return abs(a - b) <= tol


def _validate_fraction(name, value):
	if value < -EPS or value > 1 + EPS:
		raise ValueError(f"{name} must be in [0, 1], got {value}")


# Alias map for config keys -> possible carrier labels in the network.
GROUP_ALIASES = {
	"nuclear": {"nuclear"},
	"thermal": {"coal", "oil"},
	"ccgt": {"ccgt"},
	"hydro": {"hydro"},
	"phs": {"PHS", "phs"},
	"wind": {"wind", "onwind", "offwindac", "offwinddc"},
	"solar": {"solar", "pv", "solarpv"},
}

STORAGEUNIT_GROUPS = {"hydro", "phs"}


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
			for k, v in item.items():
				merged[k] = v
		return merged
	raise ValueError(f"{name} must be a mapping, got {type(value).__name__}")


class Edge:
	def __init__(self, to, rev, cap):
		self.to = to
		self.rev = rev
		self.cap = float(cap)
		self.orig_cap = float(cap)


class Dinic:
	def __init__(self, n):
		self.n = n
		self.g = [[] for _ in range(n)]

	def add_edge(self, u, v, cap):
		if cap < -EPS:
			raise ValueError(f"Negative capacity on edge {u}->{v}: {cap}")
		cap = max(0.0, float(cap))
		fwd = Edge(v, len(self.g[v]), cap)
		rev = Edge(u, len(self.g[u]), 0.0)
		self.g[u].append(fwd)
		self.g[v].append(rev)
		return len(self.g[u]) - 1

	def maxflow(self, s, t):
		flow = 0.0
		while True:
			level = [-1] * self.n
			level[s] = 0
			q = deque([s])
			while q:
				u = q.popleft()
				for e in self.g[u]:
					if e.cap > EPS and level[e.to] < 0:
						level[e.to] = level[u] + 1
						q.append(e.to)
			if level[t] < 0:
				break

			it = [0] * self.n

			def dfs(u, f):
				if u == t:
					return f
				while it[u] < len(self.g[u]):
					e = self.g[u][it[u]]
					if e.cap > EPS and level[e.to] == level[u] + 1:
						pushed = dfs(e.to, min(f, e.cap))
						if pushed > EPS:
							e.cap -= pushed
							self.g[e.to][e.rev].cap += pushed
							return pushed
					it[u] += 1
				return 0.0

			while True:
				pushed = dfs(s, float("inf"))
				if pushed <= EPS:
					break
				flow += pushed
		return flow


def _solve_bipartite_with_region_bounds(
	group_names,
	region_ids,
	group_targets,
	region_lowers,
	region_uppers,
	cell_caps,
):
	source = 0
	g_start = 1
	r_start = g_start + len(group_names)
	sink = r_start + len(region_ids)
	base_nodes = sink + 1
	super_source = base_nodes
	super_sink = base_nodes + 1

	group_node = {g: g_start + i for i, g in enumerate(group_names)}
	region_node = {r: r_start + i for i, r in enumerate(region_ids)}

	dinic = Dinic(base_nodes + 2)
	node_balance = [0.0 for _ in range(base_nodes + 2)]
	tracked_edges = {}

	def add_bounded_edge(u, v, lower, upper, tag=None):
		if upper + EPS < lower:
			raise ValueError(f"Infeasible bounds for edge {u}->{v}: [{lower}, {upper}]")
		lower = float(max(0.0, lower))
		upper = float(max(lower, upper))
		edge_idx = dinic.add_edge(u, v, upper - lower)
		node_balance[u] -= lower
		node_balance[v] += lower
		if tag is not None:
			tracked_edges[tag] = (u, edge_idx, lower)

	for g in group_names:
		add_bounded_edge(source, group_node[g], group_targets[g], group_targets[g])

	for g in group_names:
		for r in region_ids:
			cap = cell_caps.get((g, r), 0.0)
			add_bounded_edge(group_node[g], region_node[r], 0.0, cap, tag=(g, r))

	for r in region_ids:
		add_bounded_edge(region_node[r], sink, region_lowers[r], region_uppers[r])

	dinic.add_edge(sink, source, float("inf"))

	required_flow = 0.0
	for i in range(base_nodes):
		bal = node_balance[i]
		if bal > EPS:
			dinic.add_edge(super_source, i, bal)
			required_flow += bal
		elif bal < -EPS:
			dinic.add_edge(i, super_sink, -bal)

	feasible_flow = dinic.maxflow(super_source, super_sink)
	if feasible_flow + 1e-6 < required_flow:
		return None

	flows = {}
	for key, (u, idx, lower) in tracked_edges.items():
		edge = dinic.g[u][idx]
		sent = edge.orig_cap - edge.cap
		flows[key] = lower + sent

	return flows


def _bus_region_mapping(network, gadm_path, configured_regions):
	gadm = gpd.read_file(gadm_path)
	if "GADM_ID" not in gadm.columns:
		raise ValueError(f"'GADM_ID' column not found in {gadm_path}")

	gadm = gadm[gadm["GADM_ID"].isin(configured_regions)].copy()
	if gadm.empty:
		raise ValueError("No configured damage_distribution regions found in GADM shapes")

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


def _max_feasible_lower_for_region(
	group_names,
	region_ids,
	group_targets,
	region_lowers,
	region_uppers,
	cell_caps,
	target_region,
):
	base = dict(region_lowers)
	lo = 0.0
	hi = float(region_uppers[target_region])

	for _ in range(28):
		mid = 0.5 * (lo + hi)
		trial_lowers = dict(base)
		trial_lowers[target_region] = mid

		if sum(trial_lowers.values()) > sum(group_targets[g] for g in group_names) + 1e-6:
			hi = mid
			continue

		trial = _solve_bipartite_with_region_bounds(
			group_names=group_names,
			region_ids=region_ids,
			group_targets=group_targets,
			region_lowers=trial_lowers,
			region_uppers=region_uppers,
			cell_caps=cell_caps,
		)
		if trial is None:
			hi = mid
		else:
			lo = mid

	return lo


def _repair_region_shares_with_slack(
	group_names,
	region_ids,
	group_targets,
	region_shares,
	total_removed_target,
	region_uppers,
	cell_caps,
):
	shares = {r: float(region_shares[r]) for r in region_ids}
	max_iterations = max(1, len(region_ids) * len(region_ids) * 6)
	iteration = 0

	while True:
		if iteration >= max_iterations:
			raise ValueError(
				"Failed to repair regional shares: iteration limit exceeded while trimming "
				"infeasible lower bounds."
			)

		region_lowers = {
			r: max(0.0, shares[r] * total_removed_target)
			for r in region_ids
		}

		# First, trim any region that exceeds its direct upper bound.
		trimmed_share_total = 0.0
		for src in sorted(region_ids):
			upper = region_uppers[src]
			if region_lowers[src] <= upper + 1e-6:
				continue
			feasible_share = 0.0 if total_removed_target <= EPS else upper / total_removed_target
			trimmed_share_total += max(0.0, shares[src] - feasible_share)
			shares[src] = feasible_share

		if trimmed_share_total > EPS:
			while trimmed_share_total > EPS:
				best_region = None
				best_slack = -1.0
				for dst in sorted(region_ids):
					upper_share = (
						0.0
						if total_removed_target <= EPS
						else region_uppers[dst] / total_removed_target
					)
					slack = upper_share - shares[dst]
					if slack > best_slack + EPS:
						best_region = dst
						best_slack = slack

				if best_region is None or best_slack <= EPS:
					raise ValueError(
						"Infeasible generation_war_damages configuration: unable to redistribute "
						"trimmed regional share because no region has remaining slack."
					)

				delta = min(trimmed_share_total, max(0.0, best_slack))
				shares[best_region] += delta
				trimmed_share_total -= delta

		region_lowers = {
			r: max(0.0, shares[r] * total_removed_target)
			for r in region_ids
		}

		trial = _solve_bipartite_with_region_bounds(
			group_names=group_names,
			region_ids=region_ids,
			group_targets=group_targets,
			region_lowers=region_lowers,
			region_uppers=region_uppers,
			cell_caps=cell_caps,
		)
		if trial is not None:
			break

		# If still infeasible, identify regions with infeasible lower bounds under coupled carrier constraints.
		infeasible = []
		for src in sorted(region_ids):
			current_low = region_lowers[src]
			if current_low <= EPS:
				continue
			max_low = _max_feasible_lower_for_region(
				group_names=group_names,
				region_ids=region_ids,
				group_targets=group_targets,
				region_lowers=region_lowers,
				region_uppers=region_uppers,
				cell_caps=cell_caps,
				target_region=src,
			)
			if max_low + 1e-6 < current_low:
				infeasible.append((src, max_low))

		if not infeasible:
			raise ValueError(
				"Infeasible generation_war_damages configuration: bounded-flow remains infeasible "
				"after regional share repair, but no trim candidates were identified."
			)

		trimmed_share_total = 0.0
		for src, max_low in infeasible:
			new_share = 0.0 if total_removed_target <= EPS else max_low / total_removed_target
			if shares[src] > new_share + EPS:
				trimmed_share_total += shares[src] - new_share
				shares[src] = new_share

		while trimmed_share_total > EPS:
			best_region = None
			best_slack = -1.0
			for dst in sorted(region_ids):
				upper_share = (
					0.0
					if total_removed_target <= EPS
					else region_uppers[dst] / total_removed_target
				)
				slack = upper_share - shares[dst]
				if slack > best_slack + EPS:
					best_region = dst
					best_slack = slack

			if best_region is None or best_slack <= EPS:
				raise ValueError(
					"Infeasible generation_war_damages configuration: unable to redistribute "
					"trimmed regional share because no region has remaining slack."
				)

			delta = min(trimmed_share_total, max(0.0, best_slack))
			shares[best_region] += delta
			trimmed_share_total -= delta

		iteration += 1

	share_sum = sum(shares.values())
	if not _nearly_equal(share_sum, 1.0, tol=1e-6):
		raise ValueError(
			"Internal share-repair error: repaired damage_distribution shares do not sum "
			f"to 1.0 (got {share_sum:.8f})."
		)

	return shares


def main(snakemake):
	n = pypsa.Network(snakemake.input["network"])

	cfg = snakemake.config.get("generation_war_damages", {})
	if not cfg.get("enable", False):
		n.export_to_netcdf(snakemake.output["network"])
		raise SystemExit()

	damage_by_carrier = _coerce_mapping(
		"generation_war_damages.damage_percentage_by_carrier",
		cfg.get("damage_percentage_by_carrier", {}),
	)
	damage_distribution = _coerce_mapping(
		"generation_war_damages.preferred_damage_distribution",
		cfg.get("preferred_damage_distribution", {}),
	)

	if not damage_by_carrier:
		raise ValueError("generation_war_damages.damage_percentage_by_carrier is empty")
	if not damage_distribution:
		raise ValueError("generation_war_damages.preferred_damage_distribution is empty")

	region_shares = {}
	for region_id, share in damage_distribution.items():
		share = float(share)
		_validate_fraction(f"generation_war_damages.preferred_damage_distribution.{region_id}", share)
		region_shares[str(region_id)] = share

	distribution_sum = sum(region_shares.values())
	if not _nearly_equal(distribution_sum, 1.0, tol=1e-6):
		raise ValueError(
			"generation_war_damages.preferred_damage_distribution values must sum to 1.0 "
			f"(got {distribution_sum:.8f})"
		)

	generators = n.generators.copy()
	storage_units = n.storage_units.copy()
	if generators.empty and storage_units.empty:
		n.export_to_netcdf(snakemake.output["network"])
		raise SystemExit()

	if "carrier" not in generators.columns or "bus" not in generators.columns or "p_nom" not in generators.columns:
		raise ValueError("Network generators must include carrier/bus/p_nom columns")
	if "carrier" not in storage_units.columns or "bus" not in storage_units.columns or "p_nom" not in storage_units.columns:
		raise ValueError("Network storage_units must include carrier/bus/p_nom columns")

	gen_carriers_norm = generators["carrier"].map(_normalize_token)
	storage_carriers_norm = storage_units["carrier"].map(_normalize_token)

	component_tables = {
		"generator": generators,
		"storage_unit": storage_units,
	}
	component_carriers = {
		"generator": gen_carriers_norm,
		"storage_unit": storage_carriers_norm,
	}

	group_targets = {}
	group_members = {}
	group_labels = {}
	group_component = {}

	for raw_group_name, share in damage_by_carrier.items():
		share = float(share)
		_validate_fraction(
			f"generation_war_damages.damage_percentage_by_carrier.{raw_group_name}",
			share,
		)

		group_key = _normalize_token(raw_group_name)
		aliases = GROUP_ALIASES.get(group_key, {group_key})
		component = "storage_unit" if group_key in STORAGEUNIT_GROUPS else "generator"
		table = component_tables[component]
		carriers_norm = component_carriers[component]

		mask = carriers_norm.isin(aliases)
		group_assets = table.index[mask]
		total_cap = float(table.loc[group_assets, "p_nom"].clip(lower=0.0).sum())
		target_remove = share * total_cap

		if target_remove > EPS and total_cap <= EPS:
			raise ValueError(
				f"Carrier group '{raw_group_name}' has positive damage share ({share}) in "
				f"{component} domain "
				"but no installed capacity in the network"
			)

		group_targets[group_key] = target_remove
		group_members[group_key] = set(group_assets)
		group_labels[group_key] = raw_group_name
		group_component[group_key] = component

	positive_groups = [g for g, t in group_targets.items() if t > EPS]
	if not positive_groups:
		n.export_to_netcdf(snakemake.output["network"])
		raise SystemExit()

	# Ensure no asset belongs to more than one configured positive-damage group in its component domain.
	seen = {}
	for g in positive_groups:
		for gen_id in group_members[g]:
			asset_key = (group_component[g], gen_id)
			if asset_key in seen:
				comp_label = "storage_unit" if group_component[g] == "storage_unit" else "generator"
				raise ValueError(
					f"Overlapping carrier mapping detected in {comp_label}: asset "
					f"'{gen_id}' matches both '{group_labels[seen[asset_key]]}' "
					f"and '{group_labels[g]}'. Adjust aliases or config keys."
				)
			seen[asset_key] = g

	bus_to_region = _bus_region_mapping(
		network=n,
		gadm_path=snakemake.input["gadm"],
		configured_regions=set(region_shares.keys()),
	)

	gen_region = generators["bus"].map(bus_to_region)
	gen_region = gen_region.where(gen_region.notna(), None)
	storage_region = storage_units["bus"].map(bus_to_region)
	storage_region = storage_region.where(storage_region.notna(), None)
	component_regions = {
		"generator": gen_region,
		"storage_unit": storage_region,
	}

	cell_caps = defaultdict(float)
	cell_assets = defaultdict(list)

	for g in positive_groups:
		component = group_component[g]
		table = component_tables[component]
		asset_region = component_regions[component]
		for gen_id in sorted(group_members[g]):
			cap = float(max(0.0, table.at[gen_id, "p_nom"]))
			if cap <= EPS:
				continue
			region = asset_region.get(gen_id)
			if region in region_shares:
				cell_caps[(g, region)] += cap
				cell_assets[(g, region)].append((component, gen_id))

	region_ids = sorted(region_shares.keys())
	total_removed_target = sum(group_targets[g] for g in positive_groups)
	if total_removed_target <= EPS:
		n.export_to_netcdf(snakemake.output["network"])
		raise SystemExit()

	for g in positive_groups:
		regional_capacity = sum(cell_caps.get((g, r), 0.0) for r in region_ids)
		if regional_capacity + 1e-6 < group_targets[g]:
			raise ValueError(
				f"Infeasible damage target for carrier group '{group_labels[g]}': "
				f"need to remove {group_targets[g]:.6f} MW, but only {regional_capacity:.6f} MW "
				"is available in configured damage_distribution regions."
			)

	# Compute regional capacities and identify zero-capacity regions
	region_caps = {}
	region_upper_limits = {}
	for r in region_ids:
		region_cap = sum(cell_caps.get((g, r), 0.0) for g in positive_groups)
		region_caps[r] = region_cap
		region_upper_limits[r] = sum(
			min(cell_caps.get((g, r), 0.0), group_targets[g]) for g in positive_groups
		)

	if all(region_upper_limits[r] <= EPS for r in region_ids):
		raise ValueError(
			"All configured damage_distribution regions have zero installed capacity "
			"for targeted carriers. Cannot proceed."
		)

	region_shares = _repair_region_shares_with_slack(
		group_names=positive_groups,
		region_ids=region_ids,
		group_targets=group_targets,
		region_shares=region_shares,
		total_removed_target=total_removed_target,
		region_uppers=region_upper_limits,
		cell_caps=cell_caps,
	)

	region_lowers = {}
	region_uppers = {}
	for r in region_ids:
		cap = region_caps[r]
		low = max(0.0, region_shares[r] * total_removed_target)
		high = min(total_removed_target, cap)

		if low > high + 1e-6:
			raise ValueError(
				f"Infeasible regional bounds for '{r}' after share repair: "
				f"lower={low:.6f} MW, upper={high:.6f} MW."
			)

		region_lowers[r] = low
		region_uppers[r] = high

	if sum(region_lowers.values()) > total_removed_target + 1e-6:
		raise ValueError(
			"Infeasible regional bounds: sum of lower bounds exceeds total removed capacity "
			f"({sum(region_lowers.values()):.6f} > {total_removed_target:.6f})."
		)

	if sum(region_uppers.values()) + 1e-6 < total_removed_target:
		raise ValueError(
			"Infeasible regional bounds: sum of upper bounds is below total removed capacity "
			f"({sum(region_uppers.values()):.6f} < {total_removed_target:.6f})."
		)

	flows = _solve_bipartite_with_region_bounds(
		group_names=positive_groups,
		region_ids=region_ids,
		group_targets=group_targets,
		region_lowers=region_lowers,
		region_uppers=region_uppers,
		cell_caps=cell_caps,
	)

	if flows is None:
		raise ValueError(
			"Infeasible generation_war_damages configuration: cannot satisfy carrier-level "
			"damage totals and repaired regional distribution."
		)

	gen_reductions = defaultdict(float)
	storage_reductions = defaultdict(float)
	for g in positive_groups:
		for r in region_ids:
			planned = flows.get((g, r), 0.0)
			if planned <= EPS:
				continue

			remaining = planned
			for component, gen_id in cell_assets.get((g, r), []):
				table = component_tables[component]
				reductions = storage_reductions if component == "storage_unit" else gen_reductions
				avail = float(max(0.0, table.at[gen_id, "p_nom"] - reductions[gen_id]))
				if avail <= EPS:
					continue
				cut = min(avail, remaining)
				reductions[gen_id] += cut
				remaining -= cut
				if remaining <= 1e-7:
					break

			if remaining > 1e-5:
				raise ValueError(
					"Internal allocation error: unable to distribute planned cell reduction "
					f"for group '{group_labels[g]}' in region '{r}' (remaining={remaining:.6f} MW)."
				)

	for gen_id, cut in gen_reductions.items():
		new_cap = float(generators.at[gen_id, "p_nom"] - cut)
		generators.at[gen_id, "p_nom"] = max(0.0, new_cap)
	for su_id, cut in storage_reductions.items():
		new_cap = float(storage_units.at[su_id, "p_nom"] - cut)
		storage_units.at[su_id, "p_nom"] = max(0.0, new_cap)

	# Validate carrier totals.
	for g in positive_groups:
		component = group_component[g]
		reductions = storage_reductions if component == "storage_unit" else gen_reductions
		removed = sum(reductions[gen_id] for gen_id in group_members[g])
		if abs(removed - group_targets[g]) > max(1e-5, 1e-6 * group_targets[g]):
			raise ValueError(
				f"Carrier target mismatch for '{group_labels[g]}': "
				f"removed {removed:.6f} MW vs target {group_targets[g]:.6f} MW"
			)

	# Validate regional shares after repair.
	region_removed = {r: 0.0 for r in region_ids}
	for gen_id, cut in gen_reductions.items():
		region = gen_region.get(gen_id)
		if region in region_removed:
			region_removed[region] += cut
	for su_id, cut in storage_reductions.items():
		region = storage_region.get(su_id)
		if region in region_removed:
			region_removed[region] += cut

	for r in region_ids:
		share = region_removed[r] / total_removed_target
		if abs(share - region_shares[r]) > 1e-5:
			raise ValueError(
				f"Regional distribution mismatch for '{r}' after share repair: "
				f"share={share:.6f}, target={region_shares[r]:.6f}"
			)

	n.generators["p_nom"] = generators["p_nom"]
	n.storage_units["p_nom"] = storage_units["p_nom"]
	n.consistency_check()
	n.export_to_netcdf(snakemake.output["network"])


if __name__ == "__main__":
	main(snakemake)
