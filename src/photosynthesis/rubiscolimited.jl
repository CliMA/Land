###############################################################################
#
# Calculate the rubisco-limited photosynthetic rates
#
###############################################################################
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




function rubisco_limited_rate!(
            photo_set::C3ParaSet{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    @unpack Km, p_i, Vcmax, Γ_star = leaf;

    leaf.Ac = Vcmax * (p_i - Γ_star) / (p_i + Km);

    return nothing
end




function rubisco_limited_rate!(
            photo_set::C4ParaSet{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    leaf.Ac = leaf.Vcmax;

    return nothing
end




function rubisco_limited_rate!(
            photo_set::C3ParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT}
) where {FT<:AbstractFloat}
    @unpack g_lc, Km, Rd, Vcmax, Γ_star = leaf;
    @unpack p_a, p_atm = envir;

    _a = Vcmax;
    _b = Vcmax * Γ_star;
    _d = Km;
    _f = p_atm / g_lc * FT(1e-6);
    _p = p_a;

    _qa = _f;
    _qb = _f*Rd - _p - _d - _a*_f;
    _qc = _a*_p - _b - Rd*(_p + _d);
    _an = lower_quadratic(_qa, _qb, _qc);

    leaf.Ac = _an + Rd;

    return nothing
end
