#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jul-15: add method for MonoMLGrassSPAC
#     2022-Jul-15: add method for MonoMLPalmSPAC
#     2022-Jul-15: add method for MonoMLTreeSPAC
#     2022-Jul-15: rename function to plant_energy! to be more accurate (ready to add other organs other leaf)
#     2022-Jul-15: add root, trunk, branch energy budgets
#     2022-Jul-26: add leaf LMA to the denominator
#     2022-Jul-27: fix the NaN temperature issue related to grass
#     2022-Jul-27: fix the unit issue related to CP and mol to kg
#     2022-Jul-27: rescale leaf energy absorption from per groud area to per leaf area
#
#######################################################################################################################################################################################################
"""
This function has two major functionalities:
- Compute marginal energy increase in each organ
- Update the temperature in each organ when time step provided

"""
function plant_energy! end


"""

    plant_energy!(spac::MonoMLGrassSPAC{FT}) where {FT<:AbstractFloat}
    plant_energy!(spac::MonoMLPalmSPAC{FT}) where {FT<:AbstractFloat}
    plant_energy!(spac::MonoMLTreeSPAC{FT}) where {FT<:AbstractFloat}

Compute the marginal energy increase in spac, given
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` type SPAC

"""
plant_energy!(spac::MonoMLGrassSPAC{FT}) where {FT<:AbstractFloat} = (
    (; AIR, CANOPY, DIM_LAYER, DIM_ROOT, LEAVES, LEAVES_INDEX, ROOTS, ROOTS_INDEX, SOIL) = spac;

    # loop through the roots and compute the total energy out
    for _i in DIM_ROOT
        ROOTS[_i].∂e∂t  = 0;
        ROOTS[_i].∂e∂t -= flow_out(ROOTS[_i]) * CP_L_MOL(FT) * ROOTS[_i].t;
        ROOTS[_i].∂e∂t += flow_in(ROOTS[_i]) * CP_L_MOL(FT) * SOIL.LAYERS[ROOTS_INDEX[_i]].t;
    end;

    # compute the energy out from roots
    _sum_er::FT = 0;
    _sum_fr::FT = 0;
    _sum_ir::FT = 0;
    for _i in 1:DIM_ROOT
        _sum_er += flow_out(ROOTS[_i]) * CP_L_MOL(FT) * ROOTS[_i].t;
        _sum_fr += flow_out(ROOTS[_i]);
        _sum_ir += min(flow_out(ROOTS[_i]), 0);
    end;

    # if flow out the roots == 0 but _sum_er != 0, distribute the energy in roots
    if (_sum_fr == 0) && (_sum_er != 0)
        for _i in DIM_ROOT
            if flow_out(ROOTS[_i]) < 0
                ROOTS[_i].∂e∂t += _sum_er * flow_out(ROOTS[_i]) / _sum_ir;
            end;
        end;
    end;

    # loop through the leaves
    for _i in 1:DIM_LAYER
        _g_be = FT(1.4) * FT(0.135) * sqrt(AIR[LEAVES_INDEX[_i]].wind / (FT(0.72) * LEAVES[_i].WIDTH));

        LEAVES[_i].∂e∂t  = 0;
        LEAVES[_i].∂e∂t += (CANOPY.RADIATION.r_net_sw[DIM_LAYER+1-_i] + CANOPY.RADIATION.r_net_lw[DIM_LAYER+1-_i]) / (CANOPY.lai / DIM_LAYER);
        LEAVES[_i].∂e∂t -= flow_out(LEAVES[_i]) * M_H₂O(FT) * latent_heat_vapor(LEAVES[_i].t);
        LEAVES[_i].∂e∂t -= flow_out(LEAVES[_i]) * CP_L_MOL(FT) * LEAVES[_i].t;
        LEAVES[_i].∂e∂t += flow_in(LEAVES[_i]) * CP_L_MOL(FT) * LEAVES[_i].t;
        LEAVES[_i].∂e∂t -= 2 * _g_be * CP_D_MOL(FT) * (LEAVES[_i].t - AIR[LEAVES_INDEX[_i]].t);
    end;

    # compute the energy into the leaves
    _sum_el::FT = 0;
    _sum_fl::FT = 0;
    _sum_il::FT = 0;
    for _i in 1:DIM_LAYER
        _sum_el += LEAVES[_i].HS.AREA * flow_in(LEAVES[_i]) * CP_L_MOL(FT) * LEAVES[_i].t;
        _sum_fl += LEAVES[_i].HS.AREA * flow_in(LEAVES[_i]);
        _sum_il += max(LEAVES[_i].HS.AREA * flow_in(LEAVES[_i]), 0);
    end;

    # if flow into the leaves == 0 but total _sum_el != 0, distribute the energy in leaves
    if (_sum_fl == 0) && (_sum_el != 0)
        for _i in 1:DIM_LAYER
            if (flow_in(LEAVES[_i]) > 0)
                LEAVES[_i].∂e∂t -= _sum_el * (LEAVES[_i].HS.AREA * flow_in(LEAVES[_i])) / _sum_il / LEAVES[_i].HS.AREA;
            end;
        end;
    end;

    # partition the energy difference to leaves
    if (_sum_fr != 0) && (_sum_fl != 0)
        for _i in 1:DIM_LAYER
            LEAVES[_i].∂e∂t += (_sum_er - _sum_el) * (LEAVES[_i].HS.AREA * flow_in(LEAVES[_i])) / _sum_fl / LEAVES[_i].HS.AREA;
        end;
    end;

    return nothing
);

plant_energy!(spac::MonoMLPalmSPAC{FT}) where {FT<:AbstractFloat} = (
    (; AIR, CANOPY, DIM_LAYER, DIM_ROOT, LEAVES, LEAVES_INDEX, ROOTS, ROOTS_INDEX, SOIL, TRUNK) = spac;

    # loop through the roots
    TRUNK.∂e∂t = 0;
    for _i in DIM_ROOT
        ROOTS[_i].∂e∂t  = 0;
        ROOTS[_i].∂e∂t -= flow_out(ROOTS[_i]) * CP_L_MOL(FT) * ROOTS[_i].t;
        ROOTS[_i].∂e∂t += flow_in(ROOTS[_i]) * CP_L_MOL(FT) * SOIL.LAYERS[ROOTS_INDEX[_i]].t;
        TRUNK.∂e∂t     += flow_out(ROOTS[_i]) * CP_L_MOL(FT) * ROOTS[_i].t;
    end;

    # loop through the roots
    TRUNK.∂e∂t -= flow_out(TRUNK) * CP_L_MOL(FT) * TRUNK.t;

    # loop through the leaves
    for _i in 1:DIM_LAYER
        _g_be = FT(1.4) * FT(0.135) * sqrt(AIR[LEAVES_INDEX[_i]].wind / (FT(0.72) * LEAVES[_i].WIDTH));

        LEAVES[_i].∂e∂t  = 0;
        LEAVES[_i].∂e∂t += (CANOPY.RADIATION.r_net_sw[DIM_LAYER+1-_i] + CANOPY.RADIATION.r_net_lw[DIM_LAYER+1-_i]) / (CANOPY.lai / DIM_LAYER);
        LEAVES[_i].∂e∂t -= flow_out(LEAVES[_i]) * M_H₂O(FT) * latent_heat_vapor(LEAVES[_i].t);
        LEAVES[_i].∂e∂t -= flow_out(LEAVES[_i]) * CP_L_MOL(FT) * LEAVES[_i].t;
        LEAVES[_i].∂e∂t += flow_in(LEAVES[_i]) * CP_L_MOL(FT) * TRUNK.t;
        LEAVES[_i].∂e∂t -= 2 * _g_be * CP_D_MOL(FT) * (LEAVES[_i].t - AIR[LEAVES_INDEX[_i]].t);
    end;

    return nothing
);

plant_energy!(spac::MonoMLTreeSPAC{FT}) where {FT<:AbstractFloat} = (
    (; AIR, BRANCHES, CANOPY, DIM_LAYER, DIM_ROOT, LEAVES, LEAVES_INDEX, ROOTS, ROOTS_INDEX, SOIL, TRUNK) = spac;

    # loop through the roots
    TRUNK.∂e∂t = 0;
    for _i in DIM_ROOT
        ROOTS[_i].∂e∂t  = 0;
        ROOTS[_i].∂e∂t -= flow_out(ROOTS[_i]) * CP_L_MOL(FT) * ROOTS[_i].t;
        ROOTS[_i].∂e∂t += flow_in(ROOTS[_i]) * CP_L_MOL(FT) * SOIL.LAYERS[ROOTS_INDEX[_i]].t;
        TRUNK.∂e∂t     += flow_out(ROOTS[_i]) * CP_L_MOL(FT) * ROOTS[_i].t;
    end;

    # loop through the roots
    TRUNK.∂e∂t -= flow_out(TRUNK) * CP_L_MOL(FT) * TRUNK.t;

    # loop through the branches
    for _i in 1:DIM_LAYER
        BRANCHES[_i].∂e∂t  = 0;
        BRANCHES[_i].∂e∂t -= flow_out(BRANCHES[_i]) * CP_L_MOL(FT) * BRANCHES[_i].t;
        BRANCHES[_i].∂e∂t += flow_in(BRANCHES[_i]) * CP_L_MOL(FT) * TRUNK.t;
    end;

    # loop through the leaves
    for _i in 1:DIM_LAYER
        _g_be = FT(1.4) * FT(0.135) * sqrt(AIR[LEAVES_INDEX[_i]].wind / (FT(0.72) * LEAVES[_i].WIDTH));

        LEAVES[_i].∂e∂t  = 0;
        LEAVES[_i].∂e∂t += (CANOPY.RADIATION.r_net_sw[DIM_LAYER+1-_i] + CANOPY.RADIATION.r_net_lw[DIM_LAYER+1-_i]) / (CANOPY.lai / DIM_LAYER);
        LEAVES[_i].∂e∂t -= flow_out(LEAVES[_i]) * M_H₂O(FT) * latent_heat_vapor(LEAVES[_i].t);
        LEAVES[_i].∂e∂t -= flow_out(LEAVES[_i]) * CP_L_MOL(FT) * LEAVES[_i].t;
        LEAVES[_i].∂e∂t += flow_in(LEAVES[_i]) * CP_L_MOL(FT) * BRANCHES[_i].t;
        LEAVES[_i].∂e∂t -= 2 * _g_be * CP_D_MOL(FT) * (LEAVES[_i].t - AIR[LEAVES_INDEX[_i]].t);
    end;

    return nothing
);


"""

    plant_energy!(spac::MonoMLGrassSPAC{FT}, δt::FT) where {FT<:AbstractFloat}
    plant_energy!(spac::MonoMLPalmSPAC{FT}, δt::FT) where {FT<:AbstractFloat}
    plant_energy!(spac::MonoMLTreeSPAC{FT}, δt::FT) where {FT<:AbstractFloat}

Compute the marginal energy increase in spac, given
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` type SPAC
- `δt` Time step

"""
plant_energy!(spac::MonoMLGrassSPAC{FT}, δt::FT) where {FT<:AbstractFloat} = (
    (; DIM_LAYER, DIM_ROOT, LEAVES, ROOTS) = spac;

    # update the temperature for roots
    for _i in 1:DIM_ROOT
        ROOTS[_i].e += ROOTS[_i].∂e∂t * δt;
        ROOTS[_i].t  = ROOTS[_i].e / (CP_L_MOL(FT) * sum(ROOTS[_i].HS.v_storage));
    end;

    # update the temperature for leaves
    for _i in 1:DIM_LAYER
        LEAVES[_i].e += LEAVES[_i].∂e∂t * δt;
        LEAVES[_i].t  = LEAVES[_i].e / (LEAVES[_i].CP * LEAVES[_i].BIO.lma * 10 + CP_L_MOL(FT) * LEAVES[_i].HS.v_storage);
    end;

    return nothing
);

plant_energy!(spac::MonoMLPalmSPAC{FT}, δt::FT) where {FT<:AbstractFloat} = (
    (; DIM_LAYER, DIM_ROOT, LEAVES, ROOTS, TRUNK) = spac;

    # update the temperature for roots
    for _i in 1:DIM_ROOT
        ROOTS[_i].e += ROOTS[_i].∂e∂t * δt;
        ROOTS[_i].t  = ROOTS[_i].e / (CP_L_MOL(FT) * sum(ROOTS[_i].HS.v_storage));
    end;

    # update the temperature for trunk
    TRUNK.e += TRUNK.∂e∂t * δt;
    TRUNK.t  = TRUNK.e / (CP_L_MOL(FT) * sum(TRUNK.HS.v_storage));

    # update the temperature for leaves
    for _i in 1:DIM_LAYER
        LEAVES[_i].e += LEAVES[_i].∂e∂t * δt;
        LEAVES[_i].t  = LEAVES[_i].e / (LEAVES[_i].CP * LEAVES[_i].BIO.lma * 10 + CP_L_MOL(FT) * LEAVES[_i].HS.v_storage);
    end;

    return nothing
);

plant_energy!(spac::MonoMLTreeSPAC{FT}, δt::FT) where {FT<:AbstractFloat} = (
    (; BRANCHES, DIM_LAYER, DIM_ROOT, LEAVES, ROOTS, TRUNK) = spac;

    # update the temperature for roots
    for _i in 1:DIM_ROOT
        ROOTS[_i].e += ROOTS[_i].∂e∂t * δt;
        ROOTS[_i].t  = ROOTS[_i].e / (CP_L_MOL(FT) * sum(ROOTS[_i].HS.v_storage));
    end;

    # update the temperature for trunk
    TRUNK.e += TRUNK.∂e∂t * δt;
    TRUNK.t  = TRUNK.e / (CP_L_MOL(FT) * sum(TRUNK.HS.v_storage));

    # update the temperature for branches and leaves
    for _i in 1:DIM_LAYER
        BRANCHES[_i].e += BRANCHES[_i].∂e∂t * δt;
        BRANCHES[_i].t  = BRANCHES[_i].e / (CP_L_MOL(FT) * sum(BRANCHES[_i].HS.v_storage));
        LEAVES[_i].e   += LEAVES[_i].∂e∂t * δt;
        LEAVES[_i].t    = LEAVES[_i].e / (LEAVES[_i].CP * LEAVES[_i].BIO.lma * 10 + CP_L_MOL(FT) * LEAVES[_i].HS.v_storage);
    end;

    return nothing
);
