###############################################################################
#
# Calculate annual profit
#
###############################################################################
"""
    annual_simulation!(
                node::SPACSimple{FT},
                weather::DataFrame,
                output::DataFrame,
                Δt::FT = FT(1)
    ) where {FT<:AbstractFloat}

Run annual simulation for a growing season, given
- `node` [`SPACSimple`] type struct
- `weather` Weather profile in a growing season
- `output` The predefined output result
- `Δt` Time period in `[h]`
"""
function annual_simulation!(
            node::SPACSimple{FT},
            weather::DataFrame,
            output::DataFrame,
            Δt::FT = FT(1)
) where {FT<:AbstractFloat}
    # 0. unpack required values
    @unpack elevation, gaba, laba, latitude, vtoj = node;

    # 1. update the environmental constants based on the node geographycal info
    ratio            = atmospheric_pressure_ratio(elevation);
    node.envir.P_AIR = ratio * P_ATM(FT);
    node.envir.P_O₂  = node.envir.P_AIR * FT(0.209);

    # 2. calculate the growing season canopy profit
    gscp   = FT(0);
    for i in eachindex( (weather).Year )
        # 2.1 read and update the hourly data
        day  ::FT = (weather).Day[i]
        hour ::FT = (weather).Hour[i]
        _tair::FT = (weather).Tair[i]
        _dair::FT = (weather).D[i]
        p_co2::FT = (weather).CO2[i] * ratio
        r_all::FT = (weather).Solar[i]
        wind ::FT = (weather).Wind[i]
        rain ::FT = (weather).Rain[i]

        node.envir.t         = _tair + T_0(FT);
        node.envir.p_H₂O_sat = saturation_vapor_pressure( node.envir.t );
        node.envir.p_CO₂     = p_co2;
        node.envir.p_H₂O     = node.envir.p_H₂O_sat - _dair * 1000;
        node.envir.rh        = node.envir.p_H₂O / node.envir.p_H₂O_sat;
        node.envir.wind      = wind;

        # 2.2 if day time
        zenith = zenith_angle(latitude, day, hour, FT(0));
        if (r_all>0) & (zenith<=85)
            # 2.2.1 update the leaf partitioning
            big_leaf_partition!(node, zenith, r_all)
            @unpack frac_sh, frac_sl = node.container2L;

            # 2.2.2 optimize flows in each layer
            optimize_flows!(node);
            leaf_gas_exchange!(node, node.opt_f_sl, node.opt_f_sh);
            flow = node.opt_f_sl + node.opt_f_sh;
            anet = frac_sl * node.container2L.cont_sl.an + frac_sh * node.container2L.cont_sh.an;

            # 2.2.3 update drought history
            pressure_profile!(node.hs, node.p_soil, node.opt_f_sl, node.opt_f_sh, frac_sl);

            # 2.2.4 pass values to DataFrame
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

        # 2.3 if night time
        else
            # 2.3.1 calculate leaf temperature
            tlef = max(200, leaf_temperature(node, r_all, FT(0)));

            # 2.3.2 calculate the gas exchange rates
            node.ps.t = tlef;
            photosystem_temperature_dependence!(node.ps.PSM, node.envir, node.ps.t);
            node.ps.p_H₂O_sat = saturation_vapor_pressure(node.ps.t);
            anet = -node.ps.PSM.r_d;

            # 2.3.3 update temperature effects and then leaf water potential
            flow = FT(0);
            # TODO more reasonable functions need to be added
            # vc_temperature_effects!(node.hs.leaf, tlef);
            plef = end_pressure(node.hs, flow);

            # 2.3.4 pass values to DataFrame
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

            output[i, "LAI_sh"] = node.lai
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

        # pass more data to DataFrame
        output[i, "T_air" ] = _tair
        output[i, "D_air" ] = _dair
        output[i, "Wind"  ] = wind
        output[i, "Rain"  ] = rain
        output[i, "Ca"    ] = p_co2
        output[i, "SWC"   ] = node.swc
        output[i, "P_soil"] = node.p_soil
        output[i, "H_sun" ] = 90 - zenith
        output[i, "A_net" ] = anet
        output[i, "E_crit"] = node.ec

        # 2.4 update soil moisture by converting flow to Kg per hour
        flow /= KG_H_2_MOL_S(FT);
        rain *= gaba * ρ_H₂O(FT) / 1000;
        soil_moisture!(node, flow-rain, Δt);
    end

    return nothing
end
