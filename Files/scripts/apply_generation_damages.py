from collections import defaultdict, deque

import geopandas as gpd
import pypsa


EPS = 1e-9
REGIONAL_SHARE_TOLERANCE = 0.70


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
	"thermal": {"thermal", "coal", "hardcoal", "oil", "lignite", "ocgt"},
	"ccgt": {"ccgt", "combinedcyclegasturbine"},
	"hydro": {"hydro", "psp", "pshp", "ror", "runofriver"},
	"wind": {"wind", "onwind", "offwindac", "offwinddc"},
	"solar": {"solar", "pv", "solarpv"},
}


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
		"generation_war_damages.damage_distribution",
		cfg.get("damage_distribution", {}),
	)

	if not damage_by_carrier:
		raise ValueError("generation_war_damages.damage_percentage_by_carrier is empty")
	if not damage_distribution:
		raise ValueError("generation_war_damages.damage_distribution is empty")

	region_shares = {}
	for region_id, share in damage_distribution.items():
		share = float(share)
		_validate_fraction(f"generation_war_damages.damage_distribution.{region_id}", share)
		region_shares[str(region_id)] = share

	distribution_sum = sum(region_shares.values())
	if not _nearly_equal(distribution_sum, 1.0, tol=1e-6):
		raise ValueError(
			"generation_war_damages.damage_distribution values must sum to 1.0 "
			f"(got {distribution_sum:.8f})"
		)

	generators = n.generators.copy()
	if generators.empty:
		n.export_to_netcdf(snakemake.output["network"])
		raise SystemExit()

	carriers_norm = generators["carrier"].map(_normalize_token)

	group_targets = {}
	group_members = {}
	group_labels = {}

	for raw_group_name, share in damage_by_carrier.items():
		share = float(share)
		_validate_fraction(
			f"generation_war_damages.damage_percentage_by_carrier.{raw_group_name}",
			share,
		)

		group_key = _normalize_token(raw_group_name)
		aliases = GROUP_ALIASES.get(group_key, {group_key})

		mask = carriers_norm.isin(aliases)
		group_gen = generators.index[mask]
		total_cap = float(generators.loc[group_gen, "p_nom"].clip(lower=0.0).sum())
		target_remove = share * total_cap

		if target_remove > EPS and total_cap <= EPS:
			raise ValueError(
				f"Carrier group '{raw_group_name}' has positive damage share ({share}) "
				"but no installed capacity in the network"
			)

		group_targets[group_key] = target_remove
		group_members[group_key] = set(group_gen)
		group_labels[group_key] = raw_group_name

	positive_groups = [g for g, t in group_targets.items() if t > EPS]
	if not positive_groups:
		n.export_to_netcdf(snakemake.output["network"])
		raise SystemExit()

	# Ensure no generator belongs to more than one configured positive-damage group.
	seen = {}
	for g in positive_groups:
		for gen_id in group_members[g]:
			if gen_id in seen:
				raise ValueError(
					"Overlapping carrier mapping detected: generator "
					f"'{gen_id}' matches both '{group_labels[seen[gen_id]]}' "
					f"and '{group_labels[g]}'. Adjust aliases or config keys."
				)
			seen[gen_id] = g

	bus_to_region = _bus_region_mapping(
		network=n,
		gadm_path=snakemake.input["gadm"],
		configured_regions=set(region_shares.keys()),
	)

	gen_region = generators["bus"].map(bus_to_region)
	gen_region = gen_region.where(gen_region.notna(), None)

	cell_caps = defaultdict(float)
	cell_generators = defaultdict(list)

	for g in positive_groups:
		for gen_id in sorted(group_members[g]):
			cap = float(max(0.0, generators.at[gen_id, "p_nom"]))
			if cap <= EPS:
				continue
			region = gen_region.get(gen_id)
			if region in region_shares:
				cell_caps[(g, region)] += cap
				cell_generators[(g, region)].append(gen_id)

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
	for r in region_ids:
		cap = sum(cell_caps.get((g, r), 0.0) for g in positive_groups)
		region_caps[r] = cap

	# Identify zero-capacity regions and renormalize shares
	zero_cap_regions = [r for r in region_ids if region_caps[r] <= EPS]
	if zero_cap_regions:
		zero_share = sum(region_shares[r] for r in zero_cap_regions)
		remaining_share = 1.0 - zero_share

		if remaining_share <= EPS:
			raise ValueError(
				"All configured damage_distribution regions have zero installed capacity "
				"for targeted carriers. Cannot proceed."
			)

		# Renormalize remaining regions
		for r in region_ids:
			if region_caps[r] > EPS:
				region_shares[r] = region_shares[r] / remaining_share

		# Filter out zero-capacity regions
		region_ids = [r for r in region_ids if region_caps[r] > EPS]

		print(f"[generation_war_damages] Removed zero-capacity regions: {zero_cap_regions}")
		print(f"[generation_war_damages] Renormalized damage_distribution shares")

	region_lowers = {}
	region_uppers = {}
	for r in region_ids:
		cap = region_caps[r]
		low = max(0.0, (region_shares[r] - REGIONAL_SHARE_TOLERANCE) * total_removed_target)
		high = min(total_removed_target, (region_shares[r] + REGIONAL_SHARE_TOLERANCE) * total_removed_target)
		high = min(high, cap)

		if low > cap + 1e-6:
			raise ValueError(
				f"Infeasible regional lower bound for '{r}': lower={low:.6f} MW, "
				f"available={cap:.6f} MW."
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
			"damage totals and regional distribution within configured tolerance."
		)

	reductions = defaultdict(float)
	for g in positive_groups:
		for r in region_ids:
			planned = flows.get((g, r), 0.0)
			if planned <= EPS:
				continue

			remaining = planned
			for gen_id in cell_generators.get((g, r), []):
				avail = float(max(0.0, generators.at[gen_id, "p_nom"] - reductions[gen_id]))
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

	for gen_id, cut in reductions.items():
		new_cap = float(generators.at[gen_id, "p_nom"] - cut)
		generators.at[gen_id, "p_nom"] = max(0.0, new_cap)

	# Validate carrier totals.
	for g in positive_groups:
		removed = sum(reductions[gen_id] for gen_id in group_members[g])
		if abs(removed - group_targets[g]) > max(1e-5, 1e-6 * group_targets[g]):
			raise ValueError(
				f"Carrier target mismatch for '{group_labels[g]}': "
				f"removed {removed:.6f} MW vs target {group_targets[g]:.6f} MW"
			)

	# Validate regional shares with tolerance.
	region_removed = {r: 0.0 for r in region_ids}
	for gen_id, cut in reductions.items():
		region = gen_region.get(gen_id)
		if region in region_removed:
			region_removed[region] += cut

	for r in region_ids:
		share = region_removed[r] / total_removed_target
		diff = abs(share - region_shares[r])
		if diff > REGIONAL_SHARE_TOLERANCE + 1e-9:
			raise ValueError(
				f"Regional distribution mismatch for '{r}': share={share:.6f}, "
				f"target={region_shares[r]:.6f}, tolerance={REGIONAL_SHARE_TOLERANCE:.6f}"
			)

	n.generators["p_nom"] = generators["p_nom"]
	n.consistency_check()
	n.export_to_netcdf(snakemake.output["network"])


if __name__ == "__main__":
	main(snakemake)
