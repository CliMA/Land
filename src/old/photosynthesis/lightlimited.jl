###############################################################################
#
# Calculate the light-limited photosynthetic rates
#
###############################################################################
"""
    light_limited_rate!(
                photo_set::Union{C3Cytochrome{FT},C3ParaSet{FT}},
                leaf::Leaf{FT}
    ) where {FT<:AbstractFloat}
    light_limited_rate!(
                photo_set::C4ParaSet{FT},
                leaf::Leaf{FT}
    ) where {FT<:AbstractFloat}
    light_limited_rate!(
                photo_set::C3ParaSet{FT},
                leaf::Leaf{FT},
                envir::AirLayer{FT}
    ) where {FT<:AbstractFloat}

Calculate the light limited photosynthetic rate, given
- `photo_set` [`C3ParaSet`](@ref) or [`C4ParaSet`](@ref) type struct
- `leaf` [`Leaf`](@ref) type struct
- `envir` [`AirLayer`](@ref) type struct
"""
function light_limited_rate!(
            photo_set::C3Cytochrome{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    @unpack J, p_i, Γ_star = leaf;
    @unpack Eff_1, Eff_2 = photo_set;

    leaf.e2c = (p_i - Γ_star) / ( Eff_1*p_i + Eff_2*Γ_star );
    leaf.Aj  = J * leaf.e2c;

    return nothing
end
