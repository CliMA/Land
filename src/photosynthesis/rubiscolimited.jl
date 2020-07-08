###############################################################################
#
# Calculate the rubisco-limited photosynthetic rates
#
###############################################################################
"""
    rubisco_limited_rate!(photo_set::C3ParaSet{FT}, leaf::Leaf{FT}) where {FT<:AbstractFloat}
    rubisco_limited_rate!(photo_set::C4ParaSet{FT}, leaf::Leaf{FT}) where {FT<:AbstractFloat}

Calculate the RubisCO limited photosynthetic rate, given
- `photo_set` [`C3ParaSet`](@ref) or [`C4ParaSet`](@ref) type struct
- `leaf` [`Leaf`](@ref) type struct
"""
function rubisco_limited_rate!(
            photo_set::C3ParaSet{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    @unpack Km, p_i, Vcmax, Γ_star = leaf

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




"""
    rubisco_limited_rate_glc!(photo_set::C3ParaSet{FT}, leaf::Leaf{FT}, envir::AirLayer{FT}) where {FT<:AbstractFloat}

Calculate the RubisCO limited photosynthetic rate from glc, given
- `photo_set` [`C3ParaSet`](@ref) or [`C4ParaSet`](@ref) type struct
- `leaf` [`Leaf`](@ref) type struct
- `envir` [`AirLayer`](@ref) type struct
"""
function rubisco_limited_rate_glc!(
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




"""
    rubisco_limited_an_glc(photo_set::C3ParaSet{FT}, leaf::Leaf{FT}, envir::AirLayer{FT}, glcs::Array{FT,1}) where {FT<:AbstractFloat}
    rubisco_limited_an_glc(photo_set::C4ParaSet{FT}, leaf::Leaf{FT}, envir::AirLayer{FT}, glcs::Array{FT,1}) where {FT<:AbstractFloat}

Calculate the RubisCO limited photosynthetic rate from glc, given
- `photo_set` [`C3ParaSet`](@ref) or [`C4ParaSet`](@ref) type struct
- `leaf` [`Leaf`](@ref) type struct
- `envir` [`AirLayer`](@ref) type struct
- `glcs` Given array of `g_lc`
"""
function rubisco_limited_an_glc(
            photo_set::C3ParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            glcs::Array{FT,1}
) where {FT<:AbstractFloat}
    @unpack Km, Rd, Vcmax, Γ_star = leaf;
    @unpack p_a, p_atm = envir;

    _a = Vcmax;
    _b = Vcmax * Γ_star;
    _d = Km;
    _f = p_atm*FT(1e-6) ./ glcs;
    _p = p_a;

    _qa = _f;
    _qb = (Rd - _a) * _f .- (_p + _d);
    _qc = _a*_p - _b - Rd*(_p + _d);
    _an = lower_quadratic.(_qa, _qb, _qc);

    return _an
end




function rubisco_limited_an_glc(
            photo_set::C4ParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            glcs::Array{FT,1}
) where {FT<:AbstractFloat}
    _an = leaf.Vcmax

    return _an
end
