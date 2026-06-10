from pathlib import Path

import pypsa


CANDIDATE_SOLVERS = ["gurobi", "highs", "glpk"]


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
                result = n.optimize(solver_name=solver)
                if isinstance(result, tuple) and len(result) >= 2:
                    status = str(result[0]).lower()
                    condition = str(result[1]).lower()
                    if "ok" in status or "optimal" in condition:
                        return solver
                else:
                    return solver
            elif hasattr(n, "lopf"):
                ok = n.lopf(pyomo=False, solver_name=solver)
                if ok is not False:
                    return solver
            else:
                raise RuntimeError("No supported PyPSA optimization method found")
        except Exception as err:  # noqa: BLE001
            last_error = err
            continue

    raise RuntimeError(f"Optimization failed for all candidate solvers: {last_error}")


def _set_line_capacity_factor(n: pypsa.Network, factor: float) -> None:
    if n.lines.empty or "s_max_pu" not in n.lines.columns:
        return
    n.lines.loc[:, "s_max_pu"] = factor
    print(f"[base_validation] Set s_max_pu = {factor} on {len(n.lines)} lines")


def _add_load_shedding(n: pypsa.Network) -> None:
    if hasattr(n, "optimize") and hasattr(n.optimize, "add_load_shedding"):
        added = n.optimize.add_load_shedding(sign=1, marginal_cost=100_000.0, p_nom=1e6)
        if len(added):
            n.generators.loc[added, "carrier"] = "load shedding"
            if "load shedding" not in n.carriers.index:
                n.add("Carrier", "load shedding", color="#c20808")
            print(f"[base_validation] Added {len(added)} load-shedding generators")
        return

    if n.loads.empty:
        return

    if "load_shedding" not in n.carriers.index:
        n.add("Carrier", "load_shedding")

    for load_name, load_row in n.loads.iterrows():
        bus = load_row.get("bus")
        if bus not in n.buses.index:
            continue
        generator_name = f"load_shed_{load_name}"
        if generator_name in n.generators.index:
            continue
        n.add(
            "Generator",
            name=generator_name,
            bus=bus,
            carrier="load shedding",
            p_nom=1e6,
            p_nom_extendable=False,
            marginal_cost=100_000.0,
            p_max_pu=1.0,
        )


def main(snakemake):
    network_in = Path(snakemake.input["network"])
    network_out = Path(snakemake.output["network"])

    line_capacity_factor = float(snakemake.config.get("line_capacity_factor", 0.7))

    n = pypsa.Network(network_in)
    _freeze_expansion(n)
    _set_storage_cyclic(n)
    _set_line_capacity_factor(n, line_capacity_factor)
    _add_load_shedding(n)
    _solve_network(n)
    n.consistency_check()

    network_out.parent.mkdir(parents=True, exist_ok=True)
    n.export_to_netcdf(network_out)
    print(f"[base_validation] Exported solved network to {network_out}")


if __name__ == "__main__":
    main(snakemake)