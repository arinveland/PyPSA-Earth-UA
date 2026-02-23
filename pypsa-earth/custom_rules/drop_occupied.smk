rule ua_drop_occupied_network:
    input:
        network="networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
        gadm="resources/" + RDIR + "shapes/gadm_shapes.geojson"
    output:
        network="networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_pruned.nc"
    script:
        "../scripts/drop_occupied.py"

