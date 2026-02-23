import pypsa

n = pypsa.Network(snakemake.input["network"])

cfg = snakemake.config.get("outages", {})
if not cfg.get("enable", False):
    n.export_to_netcdf(snakemake.output["network"])
    raise SystemExit()

# YAML may parse numeric-looking IDs as int; make them strings to match bus IDs
drop_buses = [str(x) for x in cfg.get("buses", [])]
drop_lines = [str(x) for x in cfg.get("lines", [])]

# remove explicitly specified lines by id
if drop_lines:
    mask = n.lines.index.astype(str).isin(drop_lines)
    if mask.any():
        n.mremove("Line", n.lines.index[mask])

# remove all lines connected to specified buses
if drop_buses and len(n.lines):
    bus_mask0 = n.lines.bus0.astype(str).isin(drop_buses)
    bus_mask1 = n.lines.bus1.astype(str).isin(drop_buses)
    connected_lines = n.lines.index[bus_mask0 | bus_mask1]
    if len(connected_lines):
        n.mremove("Line", connected_lines)

# remove components on those buses, but keep buses and loads
if drop_buses:
    drop_bus_mask = n.buses.index.astype(str).isin(drop_buses)
    buses_flagged = set(n.buses.index[drop_bus_mask].astype(str))

    def drop_by_bus(component, bus_col="bus"):
        """Remove rows of `component` whose `bus_col` is in drop_buses."""
        df = getattr(n, component.lower() + "s", None)
        if df is None or bus_col not in df.columns:
            return
        mask = df[bus_col].astype(str).isin(buses_flagged)
        names = df.index[mask]
        if len(names):
            n.mremove(component, names)

    # Single-bus components: DO NOT remove Loads
    for comp in ["Generator", "StorageUnit", "Store", "ShuntImpedance"]:
        if hasattr(n, comp.lower() + "s"):
            drop_by_bus(comp, "bus")

    # Two-bus components: Links (lines already handled above)
    if hasattr(n, "links"):
        mask = (
            n.links.bus0.astype(str).isin(buses_flagged)
            | n.links.bus1.astype(str).isin(buses_flagged)
        )
        if mask.any():
            n.mremove("Link", n.links.index[mask])

    # NOTE: we intentionally DO NOT remove the buses themselves
    # and DO NOT remove n.loads – so demand remains and shows up
    # as lost load when the system cannot serve it.

n.export_to_netcdf(snakemake.output["network"])
