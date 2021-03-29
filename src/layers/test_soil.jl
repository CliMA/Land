# These functions are for testing purposes only
# They need to be tested and abstractized

# remember that the psoils are corrected by surface tension, psoil = p_25 * f_st
#=
function test_soil_from_psoil(
            node::SPACMono{FT},
            psoils::Array{FT,1},
            Δt::FT=FT(10)
) where {FT<:AbstractFloat}
    # unpack the values
    @unpack n_canopy = node;
    n_sl = node.canopy_rt.nAzi * node.canopy_rt.nIncl;

    # update the soil water contents and potential in each layer
    for i_root in 1:node.n_root
        node.swc[i_root] = max(node.mswc[i_root], soil_swc(node.plant_hs.roots[i_root].sh, psoils[i_root]));
        node.plant_hs.roots[i_root].p_ups = psoils[i_root];
    end

    e_sum = FT(0);

    # loop through the canopy layers
    for i_can in 1:node.n_canopy
        rt_layer = node.n_canopy + 1 - i_can;
        f_view = (node.can_opt.Ps[rt_layer] + node.can_opt.Ps[rt_layer+1]) / 2;

        e_i_can = FT(0);
        for i_leaf in 1:(n_sl+1)
            e_i_can += node.plant_ps[i_can].g_lw[i_leaf] *
                       (node.plant_ps[i_can].p_sat - node.envirs[i_can].p_H₂O) /
                       node.envirs[i_can].p_atm *
                       node.plant_ps[i_can].LAIx[i_leaf] *
                       node.plant_ps[i_can].LA;
        end

        # calculate the photosynthetic rates
        gas_exchange!(node.photo_set, node.plant_ps[i_can], node.envirs[i_can], GswDrive());
        gsw_control!(node.photo_set, node.plant_ps[i_can], node.envirs[i_can]);

        # use the ball-berry model here for now as the ∂A/∂E and ∂A/∂Θ functions are not yet ready
        gsw_ss = stomatal_conductance(node.stomata_model, node.plant_ps[i_can], node.envirs[i_can], FT(1));

        # assume τ = 10 minutes
        for i_leaf in 1:(n_sl+1)
            node.plant_ps[i_can].g_sw[i_leaf] += (gsw_ss[i_leaf] - node.plant_ps[i_can].g_sw[i_leaf]) / 600 * Δt;
        end

        e_sum += e_i_can;
    end

    # update root flow rates
    roots_flow!(node.plant_hs.roots,
                node.plant_hs.container_k,
                node.plant_hs.container_p,
                node.plant_hs.container_q,
                e_sum);

    return nothing
end




function test_soil_from_swc(
            node::SPACMono{FT},
            swcs::Array{FT,1},
            Δt::FT=FT(10)
) where {FT<:AbstractFloat}
    # unpack the values
    @unpack n_canopy = node;
    n_sl = node.canopy_rt.nAzi * node.canopy_rt.nIncl;

    # update the soil water contents and potential in each layer
    for i_root in 1:node.n_root
        node.swc[i_root] = max(node.mswc[i_root], swcs[i_root]);
        node.p_soil[i_root] = soil_p_25_swc(node.plant_hs.roots[i_root].sh,
                                            node.swc[i_root]) *
                              node.plant_hs.roots[i_root].f_st;
        # pass the pressure to root_hs
        node.plant_hs.roots[i_root].p_ups = node.p_soil[i_root];
    end

    e_sum = FT(0);

    # loop through the canopy layers
    for i_can in 1:node.n_canopy
        rt_layer = node.n_canopy + 1 - i_can;
        f_view = (node.can_opt.Ps[rt_layer] + node.can_opt.Ps[rt_layer+1]) / 2;

        e_i_can = FT(0);
        for i_leaf in 1:(n_sl+1)
            e_i_can += node.plant_ps[i_can].g_lw[i_leaf] *
                       (node.plant_ps[i_can].p_sat - node.envirs[i_can].p_H₂O) /
                       node.envirs[i_can].p_atm *
                       node.plant_ps[i_can].LAIx[i_leaf] *
                       node.plant_ps[i_can].LA;
        end

        # calculate the photosynthetic rates
        gas_exchange!(node.photo_set, node.plant_ps[i_can], node.envirs[i_can], GswDrive());
        gsw_control!(node.photo_set, node.plant_ps[i_can], node.envirs[i_can]);

        # use the ball-berry model here for now as the ∂A/∂E and ∂A/∂Θ functions are not yet ready
        gsw_ss = stomatal_conductance(node.stomata_model, node.plant_ps[i_can], node.envirs[i_can], FT(1));

        # assume τ = 10 minutes
        for i_leaf in 1:(n_sl+1)
            node.plant_ps[i_can].g_sw[i_leaf] += (gsw_ss[i_leaf] - node.plant_ps[i_can].g_sw[i_leaf]) / 600 * Δt;
        end

        e_sum += e_i_can;
    end

    # update root flow rates
    roots_flow!(node.plant_hs.roots,
                node.plant_hs.container_k,
                node.plant_hs.container_p,
                node.plant_hs.container_q,
                e_sum);

    return nothing
end
=#




function test_soil_from_psoil!(
            node::SPACMono{FT},
            psoils::Array{FT,1},
            e_sum::FT
) where {FT<:AbstractFloat}
    # update the soil water contents and potential in each layer
    for i_root in 1:node.n_root
        node.swc[i_root] = max(node.mswc[i_root],
                               soil_swc(node.plant_hs.roots[i_root].sh,
                                        psoils[i_root])
                              );
        node.plant_hs.roots[i_root].p_ups = psoils[i_root];
    end

    # update root flow rates
    roots_flow!(node.plant_hs, e_sum);

    return nothing
end




function test_soil_from_swc!(
            node::SPACMono{FT},
            swcs::Array{FT,1},
            e_sum::FT
) where {FT<:AbstractFloat}
    # update the soil water contents and potential in each layer
    for i_root in 1:node.n_root
        node.swc[i_root] = max(node.mswc[i_root], swcs[i_root]);
        node.p_soil[i_root] = soil_p_25_swc(node.plant_hs.roots[i_root].sh,
                                            node.swc[i_root]) *
                              node.plant_hs.roots[i_root].f_st;
        # pass the pressure to root_hs
        node.plant_hs.roots[i_root].p_ups = node.p_soil[i_root];
    end

    # update root flow rates
    roots_flow!(node.plant_hs, e_sum);

    return nothing
end
