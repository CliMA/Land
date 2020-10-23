###############################################################################
#
# Calculate the light-limited photosynthetic rates
#
###############################################################################
"""
    light_limited_rate!(
                photo_set::AbstractPhotoModelParaSet{FT},
                leaf::Leaf{FT}) where {FT<:AbstractFloat}

Calculate the Light limited photosynthetic rate, given
- `photo_set` [`C3ParaSet`](@ref) or [`C4ParaSet`](@ref) type struct
- `leaf` [`Leaf`](@ref) type struct
"""
function light_limited_rate!(
            photo_set::C3ParaSet{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    @unpack J, p_i, Γ_star = leaf;
    @unpack Eff_1, Eff_2 = photo_set;

    leaf.CO₂_per_electron = (p_i - Γ_star) / ( Eff_1*p_i + Eff_2*Γ_star );
    leaf.Aj               = J * leaf.CO₂_per_electron;

    return nothing
end




function light_limited_rate!(
            photo_set::C4ParaSet{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    leaf.CO₂_per_electron = FT(1/6);
    leaf.Aj               = leaf.J / 6;

    return nothing
end




"""
    light_limited_rate_glc!(
                photo_set::C3ParaSet{FT},
                leaf::Leaf{FT},
                envir::AirLayer{FT}) where {FT<:AbstractFloat}

Calculate the Light limited photosynthetic rate from glc, given
- `photo_set` [`C3ParaSet`](@ref) type struct
- `leaf` [`Leaf`](@ref) type struct
- `envir` [`AirLayer`](@ref) type struct
"""
function light_limited_rate_glc!(
            photo_set::C3ParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT}
) where {FT<:AbstractFloat}
    @unpack g_lc, J, Rd, Γ_star = leaf;
    @unpack p_a, p_atm = envir;
    @unpack Eff_1, Eff_2 = photo_set;

    _a = J;
    _b = J * Γ_star;
    _c = Eff_1;
    _d = Eff_2*Γ_star;
    _f = p_atm / g_lc * FT(1e-6);
    _p = p_a;

    _qa = _c * _f;
    _qb = _c*_f*Rd - _c*_p - _d - _a*_f;
    _qc = _a*_p - _b - Rd*(_c*_p + _d);
    _an = lower_quadratic(_qa, _qb, _qc);

    leaf.Aj = _an + Rd;

    return nothing
end
