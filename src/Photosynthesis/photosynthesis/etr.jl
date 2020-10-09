###############################################################################
#
# Calculate the electron transport rate
#
###############################################################################
"""
    leaf_ETR!(photo_set::C3ParaSet{FT}, leaf::Leaf{FT}) where {FT<:AbstractFloat}
    leaf_ETR!(photo_set::C4ParaSet{FT}, leaf::Leaf{FT}) where {FT<:AbstractFloat}

Update the electron transport variables in the leaf struct, given
- `photo_set` [`C3ParaSet`](@ref) or [`C4ParaSet`](@ref) type struct
- `leaf` [`Leaf`](@ref) type struct
"""
function leaf_ETR!(
            photo_set::C3ParaSet{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    @unpack APAR, maxPSII, Jmax, PSII_frac = leaf;

    _Jp = PSII_frac * maxPSII * APAR;
    _J  = min(_Jp, Jmax);

    leaf.J_pot = _Jp;
    leaf.J     = _J;

    return nothing
end




function leaf_ETR!(
            photo_set::C4ParaSet{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    @unpack APAR, maxPSII, PSII_frac = leaf

    _Jp = PSII_frac * maxPSII * APAR;

    leaf.J_pot = _Jp;
    leaf.J     = _Jp;

    return nothing
end
