function Yujie111GetAnnualProfit(
            node::SPACSimple{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            weather::DataFrame,
            maxv=100.0
) where {FT<:AbstractFloat}
    row   = length((weather).Year)
    lst   = 0
    ratio = atmospheric_pressure_ratio(node.d_alti)
    node.envir.p_atm = ratio * FT(P_ATM);
    node.envir.p_O₂  = node.envir.p_atm * FT(0.209)

    # 1. repeat each growing season with a hour interval
    a_ann = 0.0
    for i in 1:row
        if (weather).Year[i]>lst
            Yujie111UpdateSoilFromSWC(node, 1.0)
            lst = (weather).Year[i]
        end
        day   = (weather).Day[i]
        hour  = (weather).Hour[i]
        _tair = (weather).Tair[i]
        _dair = (weather).D[i]
        p_co2 = (weather).CO2[i] * ratio
        r_all = (weather).Solar[i]
        wind  = (weather).Wind[i]

        # update the leaf partitioning
        zenith = zenith_angle(node.d_lati, day, hour)
        big_leaf_partition!(node, zenith, r_all)
        @unpack frac_sh, frac_sl = node.container2L;

        node.envir.t_air = _tair + K_0;
        node.envir.p_sat = saturation_vapor_pressure( node.envir.t_air );
        node.envir.p_a   = p_co2;
        node.envir.vpd   = _dair * 1000;
        node.envir.p_H₂O = node.envir.p_sat - node.envir.vpd;
        node.envir.RH    = node.envir.p_H₂O / node.envir.p_sat;
        node.envir.wind  = wind;

        if (r_all>0) & (zenith<=85)
            Yujie111GetOptimalFs(node, photo_set, zenith, r_all)
            Yujie111GetPACGTs(node, photo_set, node.opt_f_sl, node.opt_f_sh)
            flow = node.opt_f_sl + node.opt_f_sh
            anet = frac_sl * node.container2L.cont_sl.an + frac_sh * node.container2L.cont_sh.an;
            hydraulic_p_profile!(node.hs, node.p_soil, node.opt_f_sl, node.opt_f_sh, frac_sl)
        else
            f_sl = 0.0
            f_sh = 0.0
            flow = 0.0
            tlef = max(200, leaf_temperature(node, r_all, FT(0)));
            node.ps.T = tlef;
            leaf_rd!(photo_set.ReT, node.ps);
            anet = -node.ps.Rd;
        end
        # update the soil
        flow /= FT(KG_H_2_MOL_S)
        rain  = (weather).Rain[i] * 0.001 * node.gaba * ρ_H₂O
        Yujie111UpdateSoil(node, flow-rain, 1.0)
        # add up annual profit
        a_ann += anet
    end
    a_ann *= node.laba * 3600.0 * 1E-6

    # 2. calculate the construction cost
    tmpv = 0.5 * ( maxv * log(maxv/(maxv-node.ps.Vcmax25)) + node.ps.Vcmax25 )
    tmpj = 1.75 * tmpv
    c_con = node.laba * node.c_cons *
           (1.0 + tmpv * node.c_vmax * 0.5 +
                  tmpj/1.75 * node.c_vmax * 0.5)

    # 3. calculate the profit and return
    profit = a_ann - c_con
    return profit
end
