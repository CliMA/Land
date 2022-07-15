#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jul-15: add method for MonoMLGrassSPAC
#     2022-Jul-15: add method for MonoMLPalmSPAC
#     2022-Jul-15: add method for MonoMLTreeSPAC
#     2022-Jul-15: rename function to plant_energy! to be more accurate (ready to add other organs other leaf)
#     2022-Jul-15: add root, trunk, branch energy budgets
#
#######################################################################################################################################################################################################
"""

    plant_energy!(spac::MonoMLGrassSPAC{FT}) where {FT<:AbstractFloat}
    plant_energy!(spac::MonoMLPalmSPAC{FT}) where {FT<:AbstractFloat}
    plant_energy!(spac::MonoMLTreeSPAC{FT}) where {FT<:AbstractFloat}

Compute the marginal energy increase in spac, given
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` type SPAC

"""
function plant_energy! end

plant_energy!(spac::MonoMLGrassSPAC{FT}) where {FT<:AbstractFloat} = (
    @unpack AIR, CANOPY, LEAVES, LEAVES_INDEX, N_CANOPY, N_ROOT, ROOTS, ROOTS_INDEX, SOIL = spac;

    # loop through the roots
    for _i in N_ROOT
        ROOTS[_i].∂e∂t  = 0;
        ROOTS[_i].∂e∂t -= flow_out(ROOTS[_i]) * ROOTS[_i].t;
        ROOTS[_i].∂e∂t += flow_in(ROOTS[_i]) * SOIL.LAYERS[ROOTS_INDEX[_i]].t;
    end;

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

plant_energy!(spac::MonoMLPalmSPAC{FT}) where {FT<:AbstractFloat} = (
    @unpack AIR, CANOPY, LEAVES, LEAVES_INDEX, N_CANOPY, N_ROOT, ROOTS, ROOTS_INDEX, SOIL, TRUNK = spac;

    # loop through the roots
    TRUNK.∂e∂t = 0;
    for _i in N_ROOT
        ROOTS[_i].∂e∂t  = 0;
        ROOTS[_i].∂e∂t -= flow_out(ROOTS[_i]) * ROOTS[_i].t;
        ROOTS[_i].∂e∂t += flow_in(ROOTS[_i]) * SOIL.LAYERS[ROOTS_INDEX[_i]].t;
        TRUNK.∂e∂t     += flow_out(ROOTS[_i]) * ROOTS[_i].t;
    end;

    # loop through the roots
    TRUNK.∂e∂t -= flow_out(TRUNK) * TRUNK.t;

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

plant_energy!(spac::MonoMLTreeSPAC{FT}) where {FT<:AbstractFloat} = (
    @unpack AIR, BRANCHES, CANOPY, LEAVES, LEAVES_INDEX, N_CANOPY, N_ROOT, ROOTS, ROOTS_INDEX, SOIL, TRUNK = spac;

    # loop through the roots
    TRUNK.∂e∂t = 0;
    for _i in N_ROOT
        ROOTS[_i].∂e∂t  = 0;
        ROOTS[_i].∂e∂t -= flow_out(ROOTS[_i]) * ROOTS[_i].t;
        ROOTS[_i].∂e∂t += flow_in(ROOTS[_i]) * SOIL.LAYERS[ROOTS_INDEX[_i]].t;
        TRUNK.∂e∂t     += flow_out(ROOTS[_i]) * ROOTS[_i].t;
    end;

    # loop through the roots
    TRUNK.∂e∂t -= flow_out(TRUNK) * TRUNK.t;

    # loop through the branches
    for _i in 1:N_CANOPY
        BRANCHES[_i].∂e∂t  = 0;
        BRANCHES[_i].∂e∂t -= flow_out(BRANCHES[_i]) * BRANCHES[_i].t;
        BRANCHES[_i].∂e∂t += flow_in(BRANCHES[_i]) * TRUNK.t;
    end;

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
