function simulate_growing_season!(
            node::Yujie111{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            weather::DataFrame,
            output::DataFrame
) where {FT<:AbstractFloat}
    row,col = size(weather)
    lst     = 0
    lai     = node.laba / node.gaba
    ratio   = atmospheric_pressure_ratio(node.d_alti)
    envir   = AirLayer{FT}();
    envir.p_atm = ratio * FT(P_ATM);
    envir.p_O₂  = envir.p_atm * FT(0.209)

    # 1. repeat each growing season with a hour interval
    for i in 1:row
        day   = weather[i,2 ]
        hour  = weather[i,3 ]
        t_air = weather[i,7 ]
        d_air = weather[i,9 ]
        p_co2 = weather[i,10] * ratio
        r_all = weather[i,4 ]
        wind  = weather[i,6 ]
        prec  = weather[i,5 ]
        # update the leaf partitioning
        zenith = zenith_angle(node.d_lati, day, hour)
        solh   = max(0.0, 90.0-zenith)

        envir.t_air = t_air;
        envir.p_sat = saturation_vapor_pressure(t_air + FT(K_0));
        envir.p_a   = p_co2;
        envir.vpd   = d_air * 1000;
        envir.p_H₂O = envir.p_sat - d_air * 1000;
        envir.RH    = envir.p_H₂O / envir.p_sat;
        envir.wind  = wind;

        if (r_all>0) & (zenith<=85)
            canopy = Yujie111GetLeafPartition(node, zenith, r_all)
            f_sl,f_sh = Yujie111GetOptimalFs(node, photo_set, envir, zenith, r_all)
            tmp_re = Yujie111GetPACGTs(node, photo_set, f_sl, f_sh, canopy, envir)
            flow = f_sl + f_sh
            anet = tmp_re[1,3]*canopy[1] + tmp_re[2,3]*canopy[4]
            Yujie111UpdateLegacy(node, f_sl, f_sh, canopy[1])

            output[i, "LAI_sl"] = canopy[1]*lai
            output[i, "PAR_sl"] = canopy[2]
            output[i, "RAD_sl"] = canopy[3]
            output[i, "E_sl"  ] = tmp_re[1,1]
            output[i, "P_sl"  ] = tmp_re[1,2]
            output[i, "A_sl"  ] = tmp_re[1,3]
            output[i, "C_sl"  ] = tmp_re[1,4]
            output[i, "G_sl"  ] = tmp_re[1,5]
            output[i, "T_sl"  ] = tmp_re[1,6]

            output[i, "LAI_sh"] = canopy[4]*lai
            output[i, "PAR_sh"] = canopy[5]
            output[i, "RAD_sh"] = canopy[6]
            output[i, "E_sh"  ] = tmp_re[2,1]
            output[i, "P_sh"  ] = tmp_re[2,2]
            output[i, "A_sh"  ] = tmp_re[2,3]
            output[i, "C_sh"  ] = tmp_re[2,4]
            output[i, "G_sh"  ] = tmp_re[2,5]
            output[i, "T_sh"  ] = tmp_re[2,6]
        else
            canopy = [0.0 0.0 0.0 1.0 0.0 0.0]
            f_sl = 0.0
            f_sh = 0.0
            flow = 0.0
            tlef = max(-50, Yujie111GetLeafTem(node, t_air, r_all, 0.0, false, wind));
            node.ps.T = tlef;
            leaf_rd!(photo_set.ReT, node.ps);
            anet = -node.ps.Rd;
            plef = Yujie111GetP(node, 0.0)

            output[i, "LAI_sl"] = FT(0)
            output[i, "PAR_sl"] = FT(0)
            output[i, "RAD_sl"] = FT(0)
            output[i, "E_sl"  ] = FT(0)
            output[i, "P_sl"  ] = plef
            output[i, "A_sl"  ] = anet
            output[i, "C_sl"  ] = p_co2
            output[i, "G_sl"  ] = FT(0)
            output[i, "T_sl"  ] = tlef

            output[i, "LAI_sh"] = lai
            output[i, "PAR_sh"] = FT(0)
            output[i, "RAD_sh"] = FT(0)
            output[i, "E_sh"  ] = FT(0)
            output[i, "P_sh"  ] = plef
            output[i, "A_sh"  ] = anet
            output[i, "C_sh"  ] = p_co2
            output[i, "G_sh"  ] = FT(0)
            output[i, "T_sh"  ] = tlef
        end
        output[i, "T_air" ] = t_air
        output[i, "D_air" ] = d_air
        output[i, "Wind"  ] = wind
        output[i, "Rain"  ] = prec
        output[i, "Ca"    ] = p_co2
        output[i, "SWC"   ] = node.c_curr
        output[i, "P_soil"] = node.p_soil
        output[i, "H_sun" ] = solh
        output[i, "A_net" ] = anet
        # update the soil
        rain = prec * 0.001 * node.gaba * 998.0
        Yujie111UpdateSoil(node, flow-rain, 1.0)
    end
end
