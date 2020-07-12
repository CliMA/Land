###############################################################################
#
# Calculate annual profit
#
###############################################################################
function annual_profit(
            node::SPACSimple{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            weather::DataFrame
) where {FT<:AbstractFloat}
    # 0. unpack required values
    @unpack c_cons, c_vmax, d_alti, d_lati, gaba, laba, maxv, vtoj = node;

    # 1. update the environmental constants based on the node geographycal info
    ratio            = atmospheric_pressure_ratio(d_alti);
    node.envir.p_atm = ratio * FT(P_ATM);
    node.envir.p_O₂  = node.envir.p_atm * FT(0.209);

    # 2. calculate the growing season canopy profit
    gscp   = FT(0);
    for i in eachindex( (weather).Year )
        # 2.1 read and update the hourly data
        day   = (weather).Day[i]
        hour  = (weather).Hour[i]
        _tair = (weather).Tair[i]
        _dair = (weather).D[i]
        p_co2 = (weather).CO2[i] * ratio
        r_all = (weather).Solar[i]
        wind  = (weather).Wind[i]
        rain  = (weather).Rain[i]

        node.envir.t_air = _tair + K_0;
        node.envir.p_sat = saturation_vapor_pressure( node.envir.t_air );
        node.envir.p_a   = p_co2;
        node.envir.vpd   = _dair * 1000;
        node.envir.p_H₂O = node.envir.p_sat - node.envir.vpd;
        node.envir.RH    = node.envir.p_H₂O / node.envir.p_sat;
        node.envir.wind  = wind;

        # 2.2 if day time
        zenith = zenith_angle(d_lati, day, hour);
        if (r_all>0) & (zenith<=85)
            # 2.2.1 update the leaf partitioning
            big_leaf_partition!(node, zenith, r_all)
            @unpack frac_sh, frac_sl = node.container2L;

            # 2.2.2 optimize flows in each layer
            optimize_flows!(node, photo_set);
            leaf_gas_exchange!(node, photo_set, node.opt_f_sl, node.opt_f_sh);
            flow = node.opt_f_sl + node.opt_f_sh;
            anet = frac_sl * (node.container2L).cont_sl.an + frac_sh * (node.container2L).cont_sh.an;

            # 2.2.3 update drought history
            hydraulic_p_profile!(node.hs, node.p_soil, node.opt_f_sl, node.opt_f_sh, frac_sl);

        # 2.3 if night time
        else
            node.ps.T = max(200, leaf_temperature(node, r_all, FT(0)));
            leaf_rd!(photo_set.ReT, node.ps);
            flow = FT(0);
            anet = -node.ps.Rd;
        end

        # 2.4 update soil moisture by converting flow to Kg per hour
        flow /= FT(KG_H_2_MOL_S);
        rain *= gaba * ρ_H₂O / 1000;
        soil_moisture!(node, flow-rain);

        # 2.5 add up the gscp
        gscp += anet
    end

    # 3. scale gscp to mol per year
    gscp *= laba * FT(1e-6) * 3600;

    # 4. compute the construction costs
    tmpv  = ( maxv * log(maxv/(maxv-node.ps.Vcmax25)) + node.ps.Vcmax25 ) / 2;
    cons  = laba * c_cons * (1 + tmpv * c_vmax);
    gscp -= cons;

    return gscp
end
