configfile: "Files/config.yaml"

rule ua_drop_occupied_network:
    input:
        network="Files/networks/base.nc",
        gadm="Files/shapes/gadm_shapes.geojson"
    output:
        network="Files/networks/pruned.nc"
    script:
        "Files/scripts/drop_occupied.py"

rule ua_apply_generation_damages:
    input:
        network="Files/networks/pruned.nc",
        gadm="Files/shapes/gadm_shapes.geojson"
    output:
        network="Files/networks/pruned_damages.nc"
    script:
        "Files/scripts/apply_generation_damages.py"

rule ua_add_cross_border_imports:
    input:
        network="Files/networks/pruned_damages.nc"
    output:
        network="Files/networks/xborder.nc"
    script:
        "Files/scripts/add_cross_border_imports.py"

rule ua_apply_outages:
    input:
        network="Files/networks/xborder.nc"
    output:
        network="Files/networks/xborder_outaged.nc"
    script:
        "Files/scripts/apply_outages.py"