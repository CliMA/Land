"""

    leaf_energy!(spac::MonoMLGrassSPAC{FT}) where {FT<:AbstractFloat}
    leaf_energy!(spac::MonoMLPalmSPAC{FT}) where {FT<:AbstractFloat}
    leaf_energy!(spac::MonoMLTreeSPAC{FT}) where {FT<:AbstractFloat}

Compute the marginal energy increase in spac, given
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` type SPAC

"""
function leaf_energy! end

leaf_energy!(spac::MonoMLGrassSPAC{FT}) where {FT<:AbstractFloat} = (
    @unpack AIR, CANOPY, LEAVES, LEAVES_INDEX, N_CANOPY, N_ROOT, ROOTS = spac;

    # compute the mean temperature from roots (weighted by flow rate)
    _sum_f::FT = 0;
    _sum_t::FT = 0;
    for _i in 1:N_ROOT
        _sum_f += flow_out(ROOTS[_i]);
        _sum_t += flow_out(ROOTS[_i]) * ROOTS[_i].t;
    end;
    _t_mean = _sum_t / _sum_f;

    # loop through the leaves
    for _i in 1:N_CANOPY
        _g_be = FT(1.4) * FT(0.135) * sqrt(AIR[LEAVES_INDEX[_i]].wind / (FT(0.72) * LEAVES[_i].WIDTH));

        LEAVES[_i].∂e∂t  = 0;
        LEAVES[_i].∂e∂t += CANOPY.RADIATION.r_net_sw + CANOPY.RADIATION.r_net_lw;
        LEAVES[_i].∂e∂t -= flow_out(LEAVES[_i]) * latent_heat_vapor(LEAVES[_i].t);
        LEAVES[_i].∂e∂t -= flow_out(LEAVES[_i]) * LEAVES[_i].t;
        LEAVES[_i].∂e∂t += flow_in(LEAVES[_i]) * _t_mean;
        LEAVES[_i].∂e∂t -= 2 * _g_be * CP_D_MOL(FT) * (LEAVES[_i].t - AIR[LEAVES_INDEX[_i]].t);
    end;

    return nothing
);

leaf_energy!(spac::MonoMLPalmSPAC{FT}) where {FT<:AbstractFloat} = (
    @unpack AIR, CANOPY, LEAVES, LEAVES_INDEX, N_CANOPY, TRUNK = spac;

    # loop through the leaves
    for _i in 1:N_CANOPY
        _g_be = FT(1.4) * FT(0.135) * sqrt(AIR[LEAVES_INDEX[_i]].wind / (FT(0.72) * LEAVES[_i].WIDTH));

        LEAVES[_i].∂e∂t  = 0;
        LEAVES[_i].∂e∂t += CANOPY.RADIATION.r_net_sw + CANOPY.RADIATION.r_net_lw;
        LEAVES[_i].∂e∂t -= flow_out(LEAVES[_i]) * latent_heat_vapor(LEAVES[_i].t);
        LEAVES[_i].∂e∂t -= flow_out(LEAVES[_i]) * LEAVES[_i].t;
        LEAVES[_i].∂e∂t += flow_in(LEAVES[_i]) * TRUNK.t;
        LEAVES[_i].∂e∂t -= 2 * _g_be * CP_D_MOL(FT) * (LEAVES[_i].t - AIR[LEAVES_INDEX[_i]].t);
    end;

    return nothing
);

leaf_energy!(spac::MonoMLTreeSPAC{FT}) where {FT<:AbstractFloat} = (
    @unpack AIR, BRANCHES, CANOPY, LEAVES, LEAVES_INDEX, N_CANOPY = spac;

    # loop through the leaves
    for _i in 1:N_CANOPY
        _g_be = FT(1.4) * FT(0.135) * sqrt(AIR[LEAVES_INDEX[_i]].wind / (FT(0.72) * LEAVES[_i].WIDTH));

        LEAVES[_i].∂e∂t  = 0;
        LEAVES[_i].∂e∂t += CANOPY.RADIATION.r_net_sw + CANOPY.RADIATION.r_net_lw;
        LEAVES[_i].∂e∂t -= flow_out(LEAVES[_i]) * latent_heat_vapor(LEAVES[_i].t);
        LEAVES[_i].∂e∂t -= flow_out(LEAVES[_i]) * LEAVES[_i].t;
        LEAVES[_i].∂e∂t += flow_in(LEAVES[_i]) * BRANCHES[_i].t;
        LEAVES[_i].∂e∂t -= 2 * _g_be * CP_D_MOL(FT) * (LEAVES[_i].t - AIR[LEAVES_INDEX[_i]].t);
    end;

    return nothing
);
