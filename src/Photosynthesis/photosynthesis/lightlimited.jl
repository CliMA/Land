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
            photo_set::Union{C3Cytochrome{FT},C3ParaSet{FT}},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    @unpack J, p_i, Γ_star = leaf;
    @unpack Eff_1, Eff_2 = photo_set;

    leaf.e2c = (p_i - Γ_star) / ( Eff_1*p_i + Eff_2*Γ_star );
    leaf.Aj  = J * leaf.e2c;

    return nothing
end




function light_limited_rate!(
            photo_set::C4ParaSet{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    leaf.e2c = 1 / 6;
    leaf.Aj  = leaf.J / 6;

    return nothing
end




function light_limited_rate!(
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

    leaf.Aj  = _an + Rd;
    leaf.e2c = leaf.Aj / leaf.J;

    return nothing
end
