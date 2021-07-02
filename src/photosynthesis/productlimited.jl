###############################################################################
#
# Calculate the product-limited photosynthetic rates
#
###############################################################################
"""
    product_limited_rate!(
                photo_set::C3ParaSet{FT},
                leaf::Leaf{FT}
    ) where {FT<:AbstractFloat}
    product_limited_rate!(
                photo_set::C4ParaSet{FT},
                leaf::Leaf{FT}
    ) where {FT<:AbstractFloat}
    product_limited_rate!(
                photo_set::C4ParaSet{FT},
                leaf::Leaf{FT},
                envir::AirLayer{FT}
    ) where {FT<:AbstractFloat}

Calculate the product limited photosynthetic rate, given
- `photo_set` [`C3ParaSet`](@ref) or [`C4ParaSet`](@ref) type struct
- `leaf` [`Leaf`](@ref) type struct
- `envir` [`AirLayer`](@ref) type struct
"""
function product_limited_rate!(
            photo_set::C3Cytochrome{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    @unpack p_i, Vcmax, Γ_star = leaf;
    @unpack Eff_1, Eff_2 = photo_set;

    leaf.Ap = leaf.Vcmax / 2;
    leaf.J_P680_p = leaf.Ac * (Eff_1*p_i + Eff_2*Γ_star) / (p_i - Γ_star);
    leaf.J_P700_p = leaf.J_P680_c * leaf.η;

    return nothing
end




function product_limited_rate!(
            photo_set::C3ParaSet{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    leaf.Ap = leaf.Vcmax / 2;

    return nothing
end




function product_limited_rate!(
            photo_set::C4ParaSet{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    leaf.Ap = leaf.Vpmax * leaf.p_i / (leaf.p_i + leaf.Kpep);

    return nothing
end




function product_limited_rate!(
            photo_set::C4ParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT}
) where {FT<:AbstractFloat}
    @unpack g_lc, Kpep, Rd, Vpmax = leaf;
    @unpack p_a, p_atm = envir;

    _a = Vpmax;
    _d = Kpep;
    _f = p_atm / g_lc * FT(1e-6);
    _p = p_a;

    _qa = _f;
    _qb = _f*Rd - _p - _d - _a*_f;
    _qc = _a*_p - Rd*(_p + _d);
    _an = lower_quadratic(_qa, _qb, _qc);

    leaf.Ap = _an + Rd;

    return nothing
end
