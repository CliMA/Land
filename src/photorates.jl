###############################################################################
#
# Calculate the electron transport rate
#
###############################################################################
"""
    leaf_ETR!(photo_set::C3ParaSet{FT}, leaf::Leaf{FT})
    leaf_ETR!(photo_set::C4ParaSet{FT}, leaf::Leaf{FT})

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








###############################################################################
#
# Calculate the photosynthetic rates
#
###############################################################################
"""
    rubisco_limited_rate!(photo_set::C3ParaSet{FT}, leaf::Leaf{FT})
    rubisco_limited_rate!(photo_set::C4ParaSet{FT}, leaf::Leaf{FT})

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
    light_limited_rate!(photo_set::C3ParaSet{FT}, leaf::Leaf{FT})
    light_limited_rate!(photo_set::C4ParaSet{FT}, leaf::Leaf{FT})

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
    product_limited_rate!(photo_set::C3ParaSet{FT}, leaf::Leaf{FT})
    product_limited_rate!(photo_set::C4ParaSet{FT}, leaf::Leaf{FT})

Calculate the Product limited photosynthetic rate, given
- `photo_set` [`C3ParaSet`](@ref) or [`C4ParaSet`](@ref) type struct
- `leaf` [`Leaf`](@ref) type struct
"""
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








###############################################################################
#
# Calculate the photosynthetic rates using gs
#
###############################################################################
"""
    rubisco_limited_rate_glc!(
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT})

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
    light_limited_rate_glc!(
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT})

Calculate the Light limited photosynthetic rate from glc, given
- `photo_set` [`C3ParaSet`](@ref) or [`C4ParaSet`](@ref) type struct
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




"""
    product_limited_rate_glc!(
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT})

Calculate the Product limited photosynthetic rate from glc, given
- `photo_set` [`C3ParaSet`](@ref) or [`C4ParaSet`](@ref) type struct
- `leaf` [`Leaf`](@ref) type struct
- `envir` [`AirLayer`](@ref) type struct
"""
function product_limited_rate_glc!(
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








###############################################################################
#
# Calculate the net assimilation rates using gs
#
###############################################################################
"""
    leaf_pot_ETR_APAR(leaf::Leaf{FT}, PARs::Array{FT,1})

Update the electron transport variables in the leaf struct, given
- `leaf` [`Leaf`](@ref) type struct
- `PARs` Given Array of PAR
"""
function leaf_pot_ETR_APAR(
            leaf::Leaf{FT},
            PARs::Array{FT,1}
            ) where {FT<:AbstractFloat}
    @unpack maxPSII, PSII_frac = leaf;

    _Jp = PSII_frac * maxPSII * PARs;

    return _Jp
end




"""
    rubisco_limited_an_glc(
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            glcs::Array{FT,1})

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
    _qb = (Rd - _a) * _f .- _p .- _d;
    _qc = _a*_p - _b - Rd*(_p + _d);
    _an = lower_quadratic.(_qa, _qb, _qc);

    return _an
end




"""
    light_limited_an_glc(
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            glcs::Array{FT,1})

Calculate the Light limited photosynthetic rate from glc, given
- `photo_set` [`C3ParaSet`](@ref) or [`C4ParaSet`](@ref) type struct
- `leaf` [`Leaf`](@ref) type struct
- `envir` [`AirLayer`](@ref) type struct
- `glcs` Given array of `g_lc`
"""
function light_limited_an_glc(
            photo_set::C3ParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            glcs::Array{FT,1},
            Js::Array{FT,1}
            ) where {FT<:AbstractFloat}
    @unpack Rd, Γ_star = leaf;
    @unpack p_a, p_atm = envir;
    @unpack Eff_1, Eff_2 = photo_set;

    _a = Js;
    _b = Js * Γ_star;
    _c = Eff_1;
    _d = Eff_2*Γ_star;
    _f = p_atm*FT(1e-6) ./ glcs;
    _p = p_a;

    _qa = _c * _f;
    _qb = (_c*Rd .- _a) * _f .- _c*_p .- _d;
    _qc = _a*_p .- _b .- Rd*(_c*_p + _d);
    _an = lower_quadratic.(_qa, _qb, _qc);

    return _an
end




"""
    product_limited_an_glc(
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            glcs::Array{FT,1})

Calculate the Product limited photosynthetic rate from glc, given
- `photo_set` [`C3ParaSet`](@ref) or [`C4ParaSet`](@ref) type struct
- `leaf` [`Leaf`](@ref) type struct
- `envir` [`AirLayer`](@ref) type struct
- `glcs` Given array of `g_lc`
"""
function product_limited_an_glc(
            photo_set::C4ParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            glcs::Array{FT,1}
            ) where {FT<:AbstractFloat}
    @unpack Kpep, Rd, Vpmax = leaf;
    @unpack p_a, p_atm = envir;

    _a = Vpmax;
    _d = Kpep;
    _f = p_atm*FT(1e-6) ./ glcs;
    _p = p_a;

    _qa = _f;
    _qb = (Rd-_a)*_f - _p - _d;
    _qc = _a*_p - Rd*(_p + _d);
    _an = lower_quadratic(_qa, _qb, _qc);

    return _an
end
