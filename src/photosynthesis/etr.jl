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




"""
    leaf_ETR_pot_APAR(leaf::Leaf{FT}, PARs::Array{FT,1}) where {FT<:AbstractFloat}

Update the electron transport variables in the leaf struct, given
- `leaf` [`Leaf`](@ref) type struct
- `PARs` Given Array of PAR
"""
function leaf_ETR_pot_APAR(
            leaf::Leaf{FT},
            PARs::Array{FT,1}
) where {FT<:AbstractFloat}
    @unpack maxPSII, PSII_frac = leaf;

    _Jp = PSII_frac * maxPSII * PARs;

    return _Jp
end




"""
    leaf_ETR_Jps(photo_set::C3ParaSet{FT}, leaf::Leaf{FT}, Jps::Array{FT,1}) where {FT<:AbstractFloat}
    leaf_ETR_Jps(photo_set::C4ParaSet{FT}, leaf::Leaf{FT}, Jps::Array{FT,1}) where {FT<:AbstractFloat}

Update the electron transport variables in the leaf struct, given
- `photo_set` [`C3ParaSet`](@ref) or [`C4ParaSet`](@ref) type struct
- `leaf` [`Leaf`](@ref) type struct
- `Jps` Given Array of potential ETR
"""
function leaf_ETR_Jps(
            photo_set::C3ParaSet{FT},
            leaf::Leaf{FT},
            Jps::Array{FT,1}
) where {FT<:AbstractFloat}
    @unpack Jmax = leaf;

    _J = min.(Jmax, Jps)

    return _J
end




function leaf_ETR_Jps(
            photo_set::C4ParaSet{FT},
            leaf::Leaf{FT},
            Jps::Array{FT,1}
) where {FT<:AbstractFloat}
    return Jps
end
