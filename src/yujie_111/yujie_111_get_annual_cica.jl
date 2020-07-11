function Yujie111GetAnnualCiCa(
            node::SPACSimple{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            weather::DataFrame
) where {FT<:AbstractFloat}
    row   = length((weather).Year)
    lst   = 0
    ratio = atmospheric_pressure_ratio(node.d_alti)
    node.envir.p_atm = ratio * FT(P_ATM);
    node.envir.p_O₂  = node.envir.p_atm * FT(0.209)

    # 1. repeat each growing season with a hour interval
    a_ann = 0.0
    c_ann = 0.0
    a_fav = 0.0
    c_fav = 0.0
    gpp   = 0.0
    judge = true
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





            # compute these within the Yujie111GetPACGTs function?
            flow = node.opt_f_sl + node.opt_f_sh
            a_ann += frac_sl * node.container2L.cont_sl.an + frac_sh * node.container2L.cont_sh.an;
            c_ann += (frac_sl * node.container2L.cont_sl.an * node.container2L.cont_sl.c +
                      frac_sh * node.container2L.cont_sh.an * node.container2L.cont_sh.c) / p_co2
            if judge
                a_fav += frac_sl * node.container2L.cont_sl.an + frac_sh * node.container2L.cont_sh.an;
                c_fav += (frac_sl * node.container2L.cont_sl.an * node.container2L.cont_sl.c +
                          frac_sh * node.container2L.cont_sh.an * node.container2L.cont_sh.c) / p_co2
            end
            gpp  += frac_sl * node.container2L.cont_sl.ag + frac_sh * node.container2L.cont_sh.ag;






            hydraulic_p_profile!(node.hs, node.p_soil, node.opt_f_sl, node.opt_f_sh, frac_sl)
        else
            flow = 0.0
        end
        # update the soil
        flow /= FT(KG_H_2_MOL_S)
        rain  = (weather).Rain[i] * 0.001 * node.gaba * ρ_H₂O
        Yujie111UpdateSoil(node, flow-rain, 1.0)
        if (node.p_soil>0.5)
            judge = false
        end
    end

    # 2. calculate the end of growing season PLC
    # use function in the PlantHydraulics.jl

    # 3. calculate the ci:ca and return











    # update more in this one
    # maybe merge with gs simulation function
    node.cica_all = c_ann/a_ann
    node.cica_fav = c_fav/a_fav


















    return nothing
end
