function Yujie111GetAnnualCiCa(
            node::Yujie111{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            weather
) where {FT<:AbstractFloat}
    row,col = size(weather)
    lst     = 0
    ratio   = atmospheric_pressure_ratio(node.d_alti)
    envir   = AirLayer{FT}();
    envir.p_atm = ratio * FT(P_ATM);
    envir.p_O₂  = envir.p_atm * FT(0.209)

    # 1. repeat each growing season with a hour interval
    a_ann = 0.0
    c_ann = 0.0
    a_fav = 0.0
    c_fav = 0.0
    judge = true
    for i in 1:row
        if weather[i,1]>lst
            Yujie111UpdateSoilFromSWC(node, 1.0)
            lst = weather[i,1]
        end
        day   = weather[i,2 ]
        hour  = weather[i,3 ]
        t_air = weather[i,7 ]
        d_air = weather[i,9 ]
        p_co2 = weather[i,10] * ratio
        r_all = weather[i,4 ]
        wind  = weather[i,6 ]
        # update the leaf partitioning
        zenith = zenith_angle(node.d_lati, day, hour)
        canopy = Yujie111GetLeafPartition(node, zenith, r_all)

        envir.t_air = t_air;
        envir.p_sat = saturation_vapor_pressure(t_air+FT(K_0));
        envir.p_a   = p_co2;
        envir.vpd   = d_air * 1000;
        envir.p_H₂O = envir.p_sat - d_air;
        envir.RH    = envir.p_H₂O / envir.p_sat;
        envir.wind  = wind;

        if (r_all>0) & (zenith<=85)
            f_sl,f_sh = Yujie111GetOptimalFs(node, photo_set, envir, zenith, r_all)
            tmp_re = Yujie111GetPACGTs(node, photo_set, f_sl, f_sh, canopy, envir)
            flow = f_sl + f_sh
            a_ann += tmp_re[1,3]*canopy[1] + tmp_re[2,3]*canopy[4]
            c_ann += (tmp_re[1,3]*canopy[1]*tmp_re[1,4] + tmp_re[2,3]*canopy[4]*tmp_re[2,4]) / p_co2
            if judge
                a_fav += tmp_re[1,3]*canopy[1] + tmp_re[2,3]*canopy[4]
                c_fav += (tmp_re[1,3]*canopy[1]*tmp_re[1,4] + tmp_re[2,3]*canopy[4]*tmp_re[2,4]) / p_co2
            end
            Yujie111UpdateLegacy(node, f_sl, f_sh, canopy[1])
        else
            flow = 0.0
        end
        # update the soil
        rain = weather[i,5] * 0.001 * node.gaba * 998.0
        Yujie111UpdateSoil(node, flow-rain, 1.0)
        if (node.p_soil>0.5)
            judge = false
        end
    end

    # 2. calculate the end of growing season PLC
    k_all = Yujie111GetKTree(node)

    # 3. calculate the ci:ca and return
    return [c_ann/a_ann c_fav/a_fav k_all]
end
