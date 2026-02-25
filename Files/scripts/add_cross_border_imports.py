import pypsa

n = pypsa.Network(snakemake.input["network"])

#Skip this rule if cross-border generation is not enabled
cfg = snakemake.config.get("cross_border", {})
if not cfg.get("enable", False):
    n.export_to_netcdf(snakemake.output["network"])
    raise SystemExit()

carrier = cfg.get("carrier", "cross-border")
buses = cfg.get("buses", [])
ntc_gw = float(cfg.get("ntc_gw", 0.0))
mc = float(cfg.get("marginal_cost", 0.0))
p_nom = float(cfg.get("p_nom_mw", 1e6))

#Throw error if cross-border capacity is negative
if ntc_gw < 0:
    raise ValueError("cross_border.ntc_gw must be >= 0")

#Skip this rule if cross-border capacity is 0
if ntc_gw == 0:
    # no imports scenario: don't add generators or a constraint
    n.export_to_netcdf(snakemake.output["network"])
    raise SystemExit()

#Throw error if no buses are assigned to cross-border generation
if not buses:
    raise ValueError("cross_border.buses must be a non-empty list")

# Throw error if cross-border buses are undefined
missing = [b for b in buses if b not in n.buses.index]
if missing:
    raise ValueError(f"cross_border.buses not found in network buses (first 20): {missing[:20]}")

# Ensure carrier exists 
if carrier not in n.carriers.index:
    n.add("Carrier", carrier)

# Add one import generator per interconnection bus
for b in buses:
    name = f"import_{b}"
    if name in n.generators.index:
        continue
    n.add(
        "Generator",
        name=name,
        bus=b,
        carrier=carrier,
        p_nom=p_nom,
        p_nom_extendable=False,
        marginal_cost=mc,
        p_min_pu=0.0,  # import-only
        p_max_pu=1.0,
    )
n.export_to_netcdf(snakemake.output["network"])