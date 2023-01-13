#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-01: migrate function from older version
#     2022-Jul-01: rename the function from gsw_control! to limit_stomatal_conductance!
#     2022-Jul-01: add method for Leaf
#     2022-Jul-01: add method for Leaves1D
#     2022-Jul-01: add method for Leaves2D
#     2022-Jul-11: deflate documentations
#
#######################################################################################################################################################################################################
"""

    limit_stomatal_conductance!(leaf::Leaf{FT}) where {FT<:AbstractFloat}
    limit_stomatal_conductance!(leaves::Leaves1D{FT}) where {FT<:AbstractFloat}
    limit_stomatal_conductance!(leaves::Leaves2D{FT}) where {FT<:AbstractFloat}

Limit stomatal conductance for H₂O for
- `leaf` `Leaf` type struct
- `leaves` `Leaves1D` type struct

"""
function limit_stomatal_conductance! end

limit_stomatal_conductance!(leaf::Leaf{FT}) where {FT<:AbstractFloat} = (
    (; G_LIMITS) = leaf;

    _ratio = relative_diffusive_coefficient(leaf.t);
    _g_min = G_LIMITS[1] * _ratio;
    _g_max = G_LIMITS[2] * _ratio;

    # if gsw is lower than the limits
    if leaf.g_H₂O_s < _g_min
        leaf.g_H₂O_s = _g_min
    end;

    # if gsw is higher than the limits
    if leaf.g_H₂O_s > _g_max
        leaf.g_H₂O_s = _g_max
    end;

    return nothing
);

limit_stomatal_conductance!(leaves::Leaves1D{FT}) where {FT<:AbstractFloat} = (
    (; G_LIMITS) = leaves;

    _ratio_sunlit = relative_diffusive_coefficient(leaves.t[1]);
    _ratio_shaded = relative_diffusive_coefficient(leaves.t[2]);
    _g_min_sunlit = G_LIMITS[1] * _ratio_sunlit;
    _g_max_sunlit = G_LIMITS[2] * _ratio_sunlit;
    _g_min_shaded = G_LIMITS[1] * _ratio_shaded;
    _g_max_shaded = G_LIMITS[2] * _ratio_shaded;

    # for sunlit leaves
    if leaves.g_H₂O_s[1] < _g_min_sunlit
        leaves.g_H₂O_s[1] = _g_min_sunlit
    end;
    if leaves.g_H₂O_s[1] > _g_max_sunlit
        leaves.g_H₂O_s[1] = _g_max_sunlit
    end;

    # for shaded leaves
    if leaves.g_H₂O_s[2] < _g_min_shaded
        leaves.g_H₂O_s[2] = _g_min_shaded
    end;
    if leaves.g_H₂O_s[2] > _g_max_shaded
        leaves.g_H₂O_s[2] = _g_max_shaded
    end;

    return nothing
);

limit_stomatal_conductance!(leaves::Leaves2D{FT}) where {FT<:AbstractFloat} = (
    (; G_LIMITS) = leaves;

    _ratio = relative_diffusive_coefficient(leaves.t);
    _g_min = G_LIMITS[1] * _ratio;
    _g_max = G_LIMITS[2] * _ratio;

    # for sunlit leaves
    for _i in eachindex(leaves.g_H₂O_s_sunlit)
        if leaves.g_H₂O_s_sunlit[_i] < _g_min
            leaves.g_H₂O_s_sunlit[_i] = _g_min
        end;
        if leaves.g_H₂O_s_sunlit[_i] > _g_max
            leaves.g_H₂O_s_sunlit[_i] = _g_max
        end;
    end;

    # for shaded leaves
    if leaves.g_H₂O_s_shaded < _g_min
        leaves.g_H₂O_s_shaded = _g_min
    end;
    if leaves.g_H₂O_s_shaded > _g_max
        leaves.g_H₂O_s_shaded = _g_max
    end;

    return nothing
);
