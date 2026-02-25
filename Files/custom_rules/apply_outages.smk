rule ua_apply_outages:
    input:
        network="networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_xborder.nc"
    output:
        network="networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_xborder_outaged.nc"
    script:
        "../scripts/apply_outages.py"
