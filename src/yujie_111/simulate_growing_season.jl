function simulate_growing_season!(
            node::SPACSimple{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            weather::DataFrame,
            output::DataFrame
) where {FT<:AbstractFloat}
    row,col = size(weather)
    lst     = 0
    lai     = node.laba / node.gaba
    ratio   = atmospheric_pressure_ratio(node.d_alti)
    node.envir.p_atm = ratio * FT(P_ATM);
    node.envir.p_O₂  = node.envir.p_atm * FT(0.209)

    # 1. repeat each growing season with a hour interval
    for i in 1:row
        day   = weather[i,2 ]
        hour  = weather[i,3 ]
        _tair = weather[i,7 ]
        _dair = weather[i,9 ]
        p_co2 = weather[i,10] * ratio
        r_all = weather[i,4 ]
        wind  = weather[i,6 ]
        prec  = weather[i,5 ]
        # update the leaf partitioning
        zenith = zenith_angle(node.d_lati, day, hour)
        solh   = max(0.0, 90.0-zenith)

        node.envir.t_air = _tair + K_0;
        node.envir.p_sat = saturation_vapor_pressure( node.envir.t_air );
        node.envir.p_a   = p_co2;
        node.envir.vpd   = _dair * 1000;
        node.envir.p_H₂O = node.envir.p_sat - node.envir.vpd;
        node.envir.RH    = node.envir.p_H₂O / node.envir.p_sat;
        node.envir.wind  = wind;

        if (r_all>0) & (zenith<=85)
            big_leaf_partition!(node, zenith, r_all)
            @unpack frac_sh, frac_sl = (node.container2L);

            Yujie111GetOptimalFs(node, photo_set, zenith, r_all)
            Yujie111GetPACGTs(node, photo_set, node.opt_f_sl, node.opt_f_sh)
            flow = node.opt_f_sl + node.opt_f_sh
            anet = frac_sl * (node.container2L).cont_sl.an + frac_sh * (node.container2L).cont_sh.an;
            hydraulic_p_profile!(node.hs, node.p_soil, node.opt_f_sl, node.opt_f_sh, frac_sl)

            output[i, "LAI_sl"] = (node.container2L).lai_sl
            output[i, "PAR_sl"] = (node.container2L).par_sl
            output[i, "RAD_sl"] = (node.container2L).rad_sl
            output[i, "E_sl"  ] = (node.container2L).cont_sl.e
            output[i, "P_sl"  ] = (node.container2L).cont_sl.p
            output[i, "An_sl" ] = (node.container2L).cont_sl.an
            output[i, "Ag_sl" ] = (node.container2L).cont_sl.ag
            output[i, "C_sl"  ] = (node.container2L).cont_sl.c
            output[i, "G_sl"  ] = (node.container2L).cont_sl.gh
            output[i, "T_sl"  ] = (node.container2L).cont_sl.t

            output[i, "LAI_sh"] = (node.container2L).lai_sh
            output[i, "PAR_sh"] = (node.container2L).par_sh
            output[i, "RAD_sh"] = (node.container2L).rad_sh
            output[i, "E_sh"  ] = (node.container2L).cont_sh.e
            output[i, "P_sh"  ] = (node.container2L).cont_sh.p
            output[i, "An_sh" ] = (node.container2L).cont_sh.an
            output[i, "Ag_sh" ] = (node.container2L).cont_sh.ag
            output[i, "C_sh"  ] = (node.container2L).cont_sh.c
            output[i, "G_sh"  ] = (node.container2L).cont_sh.gh
            output[i, "T_sh"  ] = (node.container2L).cont_sh.t
        else
            f_sl = 0.0
            f_sh = 0.0
            flow = 0.0
            tlef = max(200, leaf_temperature(node, r_all, flow));
            node.ps.T = tlef;
            leaf_rd!(photo_set.ReT, node.ps);
            anet = -node.ps.Rd;
            plef = xylem_p_from_flow(node.hs, flow)

            output[i, "LAI_sl"] = FT(0)
            output[i, "PAR_sl"] = FT(0)
            output[i, "RAD_sl"] = FT(0)
            output[i, "E_sl"  ] = FT(0)
            output[i, "P_sl"  ] = plef
            output[i, "An_sl" ] = anet
            output[i, "Ag_sl" ] = FT(0)
            output[i, "C_sl"  ] = p_co2
            output[i, "G_sl"  ] = FT(0)
            output[i, "T_sl"  ] = tlef

            output[i, "LAI_sh"] = lai
            output[i, "PAR_sh"] = FT(0)
            output[i, "RAD_sh"] = FT(0)
            output[i, "E_sh"  ] = FT(0)
            output[i, "P_sh"  ] = plef
            output[i, "An_sh" ] = anet
            output[i, "Ag_sh" ] = FT(0)
            output[i, "C_sh"  ] = p_co2
            output[i, "G_sh"  ] = FT(0)
            output[i, "T_sh"  ] = tlef
        end
        output[i, "T_air" ] = _tair
        output[i, "D_air" ] = _dair
        output[i, "Wind"  ] = wind
        output[i, "Rain"  ] = prec
        output[i, "Ca"    ] = p_co2
        output[i, "SWC"   ] = node.c_curr
        output[i, "P_soil"] = node.p_soil
        output[i, "H_sun" ] = solh
        output[i, "A_net" ] = anet

        # convert flow into kg update the soil
        flow /= FT(KG_H_2_MOL_S)
        rain  = prec * 0.001 * node.gaba * ρ_H₂O

        Yujie111UpdateSoil(node, flow-rain, 1.0)
    end
end
