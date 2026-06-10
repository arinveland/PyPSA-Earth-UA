configfile: "Files/config.yaml"

scenario_cfg = config.get("scenario_analysis", {})
scenario_rdir = str(scenario_cfg.get("rdir_name", "")).strip()
scenario_results_dir = "Files/results" if not scenario_rdir else f"Files/results/{scenario_rdir}"

rule ua_drop_occupied_network:
    input:
        network="Files/networks/base.nc",
        gadm="Files/shapes/gadm_shapes.geojson"
    output:
        network="Files/networks/pruned.nc"
    script:
        "Files/scripts/drop_occupied.py"

rule ua_base_validation:
    input:
        network="Files/networks/pruned.nc"
    output:
        network="Files/networks/base_validation.nc"
    script:
        "Files/scripts/base_validation.py"

rule ua_adjust_loads:
    input:
        network="Files/networks/pruned.nc",
        gadm="Files/shapes/gadm_shapes.geojson"
    output:
        network="Files/networks/pruned_loads.nc"
    script:
        "Files/scripts/adjust_loads.py"

rule ua_add_cross_border_imports:
    input:
        network="Files/networks/pruned_loads.nc"
    output:
        network="Files/networks/xborder.nc"
    script:
        "Files/scripts/add_cross_border_imports.py"

rule ua_scenario_analysis:
    input:
        network="Files/networks/xborder.nc",
        gadm="Files/shapes/gadm_shapes.geojson"
    output:
        network=f"{scenario_results_dir}/scenario_analysis.nc"
    script:
        "Files/scripts/scenario_analysis.py"