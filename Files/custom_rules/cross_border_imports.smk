rule ua_add_cross_border_imports:
    input:
        network="networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_pruned.nc"
    output:
        network="networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_xborder.nc"
    script:
        "../scripts/add_cross_border_imports.py"

