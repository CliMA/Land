###############################################################################
#
# Calculate the rubisco-limited photosynthetic rates
#
###############################################################################
#=
"""
    rubisco_limited_rate!(
                photo_set::Union{C3Cytochrome{FT},C3ParaSet{FT}},
                leaf::Leaf{FT}
    ) where {FT<:AbstractFloat}
    rubisco_limited_rate!(
                photo_set::C4ParaSet{FT},
                leaf::Leaf{FT}
    ) where {FT<:AbstractFloat}
    rubisco_limited_rate!(
                photo_set::C3ParaSet{FT},
                leaf::Leaf{FT},
                envir::AirLayer{FT}
    ) where {FT<:AbstractFloat}

Calculate the RubisCO limited photosynthetic rate, given
- `photo_set` [`C3ParaSet`](@ref) or [`C4ParaSet`](@ref) type struct
- `leaf` [`Leaf`](@ref) type struct
- `envir` [`AirLayer`](@ref) type struct
"""
function rubisco_limited_rate!(
            photo_set::C3Cytochrome{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    @unpack Km, p_i, Vcmax, Γ_star = leaf;
    @unpack Eff_1, Eff_2 = photo_set;

    leaf.Ac = Vcmax * (p_i - Γ_star) / (p_i + Km);
    leaf.J_P680_c = leaf.Ac * (Eff_1*p_i + Eff_2*Γ_star) / (p_i - Γ_star);
    leaf.J_P700_c = leaf.J_P680_c * leaf.η;

    return nothing
end
=#
