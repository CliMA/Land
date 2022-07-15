
function leaf_energy! end

leaf_energy!(spac::MonoMLTreeSPAC{FT}) where {FT<:AbstractFloat} = (
    @unpack AIR, BRANCHES, CANOPY, LEAVES, LEAVES_INDEX, N_CANOPY = spac;

    # loop through the leaves
    for _i in 1:N_CANOPY
        _g_be = FT(1.4) * FT(0.135) * sqrt(AIR[LEAVES_INDEX[_i]].wind / (FT(0.72) * LEAVES[_i].WIDTH));

        LEAVES[_i].∂e∂t  = 0;
        LEAVES[_i].∂e∂t += CANOPY.RADIATION.r_net_sw + CANOPY.RADIATION.r_net_lw;
        LEAVES[_i].∂e∂t -= leaf_sink(LEAVES[_i]) * latent_heat_vapor(LEAVES[_i].t);
        LEAVES[_i].∂e∂t += leaf_source(LEAVES[_i]) * BRANCHES[_i].t;
        LEAVES[_i].∂e∂t -= 2 * _g_be * CP_D_MOL(FT) * (LEAVES[_i].t - AIR[LEAVES_INDEX[_i]].t);
    end;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-15: add function to read net water entering the leaf
#
#######################################################################################################################################################################################################
"""

    leaf_sink(lf::Union{Leaf{FT}, Leaves2D{FT}}) where {FT<:AbstractFloat}

Return the net flow that enter the leaf, given
- `lf` `Leaf` or `Leaves2D` type leaf

"""
function leaf_sink end

leaf_sink(lf::Union{Leaf{FT}, Leaves2D{FT}}) where {FT<:AbstractFloat} = leaf_sink(lf.HS.FLOW);

leaf_sink(mode::SteadyStateFlow{FT}) where {FT<:AbstractFloat} = FT(0);

leaf_sink(mode::NonSteadyStateFlow{FT}) where {FT<:AbstractFloat} = mode.f_in - mode.f_out;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-15: add function to read water exiting the leaf
#
#######################################################################################################################################################################################################
"""

    leaf_source(lf::Union{Leaf{FT}, Leaves2D{FT}}) where {FT<:AbstractFloat}

Return the net flow that escape from the leaf, given
- `lf` `Leaf` or `Leaves2D` type leaf

"""
function leaf_source end

leaf_source(lf::Union{Leaf{FT}, Leaves2D{FT}}) where {FT<:AbstractFloat} = leaf_source(lf.HS.FLOW);

leaf_source(mode::SteadyStateFlow{FT}) where {FT<:AbstractFloat} = mode.flow;

leaf_source(mode::NonSteadyStateFlow{FT}) where {FT<:AbstractFloat} = mode.f_out;
