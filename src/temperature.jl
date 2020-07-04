###############################################################################
#
# Arrhenius corrections
#
###############################################################################
"""
    arrhenius_correction(td_set::AbstractTDParameterSet{FT}, T::FT)

A correction factor based on arrhenius's fitting procedure, given
- `td_set` [`ArrheniusTD`](@ref) or [`ArrheniusPeakTD`](@ref) type struct
- `T` Leaf temperature in `[K]`

The equation used for [`ArrheniusTD`](@ref) is
```math
corr = \\exp \\left( \\dfrac{ΔHa}{R T_0} - \\dfrac{ΔHa}{R T_1} \\right)
```

The equations used for [`ArrheniusPeakTD`](@ref) are
```math
corr = \\exp \\left( \\dfrac{ΔHa}{R T_0} - \\dfrac{ΔHa}{R T_1} \\right)
       \\cdot
       \\dfrac{ 1 + \\exp \\left( \\dfrac{S_v T_0 - H_d}{R T_0} \\right) }
              { 1 + \\exp \\left( \\dfrac{S_v T_1 - H_d}{R T_1} \\right) }
```
"""
function arrhenius_correction(
            td_set::ArrheniusTD{FT},
            T::FT
            ) where {FT<:AbstractFloat}
    return exp( td_set.ΔHa_to_RT25 - td_set.ΔHa_to_R/T )
end

function arrhenius_correction(
            td_set::ArrheniusPeakTD{FT},
            T::FT
            ) where {FT<:AbstractFloat}
    @unpack C, ΔHa_to_RT25, ΔHd_to_R, ΔSv_to_R = td_set

    # _f_a: activation correction, C/_f_b: de-activation correction
    _f_a::FT = exp( ΔHa_to_RT25 * (1 - FT(K_25)/T) );
    _f_b::FT = 1 + exp(ΔSv_to_R - ΔHd_to_R/T);

    return C /_f_b * _f_a
end








###############################################################################
#
# Temperature dependency of the photosynthetic parameters
#
###############################################################################
"""
    photo_TD_from_set(td_set::ArrheniusTD{FT}, T::FT)

Make temperature correction from parameter set, given
- `td_set` [`ArrheniusTD`](@ref) type parameter set, which has a `VAL_25` field
- `T` Leaf temperature

Useful for Kc, Ko, Kpep, and ``Γ^{*}``.
"""
function photo_TD_from_set(
            td_set::ArrheniusTD{FT},
            T::FT
            ) where {FT<:AbstractFloat}
    return td_set.VAL_25 * arrhenius_correction(td_set, T)
end




"""
    photo_TD_from_val(td_set::AbstractTDParameterSet, val::FT, T::FT)

Make temperature correction from a given value, given
- `td_set` [`ArrheniusTD`](@ref) or [`ArrheniusPeakTD`](@ref) type struct
- `val` Uncorrected value at 298.15 K
- `T` Leaf temperature

Useful for Vcmax, Vomac, Vpmax, Jmax, and Respiration.
"""
function photo_TD_from_val(
            td_set::AbstractTDParameterSet{FT},
            val::FT,
            T::FT
            ) where {FT<:AbstractFloat}
    return val * arrhenius_correction(td_set, T)
end








###############################################################################
#
# Functions to update the TD individually
#
###############################################################################
"""
    leaf_jmax!(td_set::AbstractTDParameterSet, leaf::Leaf{FT})

Update maximal electron transport rate at leaf temperature, given
- `td_set` [`AbstractTDParameterSet`](@ref) type TD parameter set
- `leaf` [`Leaf`](@ref) type struct
"""
function leaf_jmax!(
            td_set::AbstractTDParameterSet{FT},
            leaf::Leaf{FT}
            ) where {FT<:AbstractFloat}
    leaf.Jmax = photo_TD_from_val(td_set, leaf.Jmax25, leaf.T);

    return nothing
end




"""
    leaf_kc!(td_set::ArrheniusTD{FT}, leaf::Leaf{FT})

Update Kc at leaf temperature, given
- `td_set` [`ArrheniusTD`](@ref) type TD parameter set
- `leaf` [`Leaf`](@ref) type struct
"""
function leaf_kc!(
            td_set::ArrheniusTD{FT},
            leaf::Leaf{FT}
            ) where {FT<:AbstractFloat}
    leaf.Kc = photo_TD_from_set(td_set, leaf.T);

    return nothing
end




"""
    leaf_km!(photo_set::C3ParaSet{FT}, leaf::Leaf{FT}, envir::AirLayer{FT})

Update Ko at leaf temperature, given
- `photo_set` [`C3ParaSet`](@ref) type photosynthesis parameter set
- `leaf` [`Leaf`](@ref) type struct
- `envir` [`AirLayer`](@ref) type struct
"""
function leaf_km!(
            photo_set::C3ParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT}
            ) where {FT<:AbstractFloat}
    leaf.Km = leaf.Kc * (1 + envir.p_O₂/leaf.Ko);

    return nothing
end




"""
    leaf_ko!(td_set::ArrheniusTD{FT}, leaf::Leaf{FT})

Update Ko at leaf temperature, given
- `td_set` [`ArrheniusTD`](@ref) type TD parameter set
- `leaf` [`Leaf`](@ref) type struct
"""
function leaf_ko!(
            td_set::ArrheniusTD{FT},
            leaf::Leaf{FT}
            ) where {FT<:AbstractFloat}
    leaf.Ko = photo_TD_from_set(td_set, leaf.T);

    return nothing
end




"""
    leaf_kpep!(td_set::ArrheniusTD{FT}, leaf::Leaf{FT})

Update Kpep at leaf temperature, given
- `td_set` [`ArrheniusTD`](@ref) type TD parameter set
- `leaf` [`Leaf`](@ref) type struct
"""
function leaf_kpep!(
            td_set::ArrheniusTD{FT},
            leaf::Leaf{FT}
            ) where {FT<:AbstractFloat}
    leaf.Kpep = photo_TD_from_set(td_set, leaf.T);

    return nothing
end




"""
    leaf_rd!(td_set::AbstractTDParameterSet, leaf::Leaf)

Update leaf dark respiration rate at leaf temperature, given
- `td_set` [`AbstractTDParameterSet`](@ref) type TD parameter set
- `leaf` [`Leaf`](@ref) type struct
"""
function leaf_rd!(
            td_set::AbstractTDParameterSet{FT},
            leaf::Leaf{FT}
            ) where {FT<:AbstractFloat}
    leaf.Rd = photo_TD_from_val(td_set, leaf.Rd25, leaf.T);

    return nothing
end




"""
    leaf_vcmax!(td_set::AbstractTDParameterSet, leaf::Leaf)

Update leaf maximal carboxylation rate at leaf temperature, given
- `td_set` [`AbstractTDParameterSet`](@ref) type TD parameter set
- `leaf` [`Leaf`](@ref) type struct
"""
function leaf_vcmax!(
            td_set::AbstractTDParameterSet{FT},
            leaf::Leaf{FT}
            ) where {FT<:AbstractFloat}
    leaf.Vcmax = photo_TD_from_val(td_set, leaf.Vcmax25, leaf.T);

    return nothing
end




"""
    leaf_vpmax!(td_set::AbstractTDParameterSet, leaf::Leaf)

Update leaf maximal PEP carboxylation rate at leaf temperature, given
- `td_set` [`AbstractTDParameterSet`](@ref) type TD parameter set
- `leaf` [`Leaf`](@ref) type struct
"""
function leaf_vpmax!(
            td_set::AbstractTDParameterSet{FT},
            leaf::Leaf{FT}
            ) where {FT<:AbstractFloat}
    leaf.Vpmax = photo_TD_from_val(td_set, leaf.Vpmax25, leaf.T);

    return nothing
end




"""
    leaf_Γstar!(td_set::ArrheniusTD{FT}, leaf::Leaf{FT})

Update ``Γ^{*}`` at leaf temperature, given
- `td_set` [`ArrheniusTD`](@ref) type TD parameter set
- `leaf` [`Leaf`](@ref) type struct
"""
function leaf_Γstar!(
            td_set::ArrheniusTD{FT},
            leaf::Leaf{FT}
            ) where {FT<:AbstractFloat}
    leaf.Γ_star = photo_TD_from_set(td_set, leaf.T);

    return nothing
end








###############################################################################
#
# Calculate photosynthesis using Leaf
#
###############################################################################
"""
    leaf_temperature_dependence!(
            photo_set::AbstractPhotoModelParaSet,
            leaf::AbstractLeaf,
            envir::AirLayer{FT})
    leaf_temperature_dependence!(
            photo_set::AbstractPhotoModelParaSet,
            leaf::AbstractLeaf,
            envir::AirLayer{FT},
            T::FT)

Update the temperature dependent photosynthesis only, given
- `photo_set` [`AbstractPhotoModelParaSet`](@ref) type parameter set
- `leaf` [`Leaf`](@ref) type struct
- `envir` [`AirLayer`](@ref) type struct
- `T` Given leaf temperature
"""
function leaf_temperature_dependence!(
            photo_set::C3ParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT}
            ) where {FT<:AbstractFloat}
    leaf.p_sat = saturation_vapor_pressure(leaf.T);
    leaf_rd!(photo_set.ReT, leaf);
    leaf_vcmax!(photo_set.VcT, leaf);
    leaf_jmax!(photo_set.JT , leaf);
    leaf_kc!(photo_set.KcT, leaf);
    leaf_ko!(photo_set.KoT, leaf);
    leaf_km!(photo_set, leaf, envir);
    leaf_Γstar!(photo_set.ΓsT, leaf);

    return nothing
end

function leaf_temperature_dependence!(
            photo_set::C4ParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT}
            ) where {FT<:AbstractFloat}
    leaf.p_sat = saturation_vapor_pressure(leaf.T);
    leaf_rd!(photo_set.ReT, leaf);
    leaf_vcmax!(photo_set.VcT, leaf);
    leaf_vpmax!(photo_set.VpT, leaf);
    leaf_kpep!(photo_set.KpT, leaf);

    return nothing
end

function leaf_temperature_dependence!(
            photo_set::C4ParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT},
            T::FT
            ) where {FT<:AbstractFloat}
    leaf.T = T;
    leaf_temperature_dependence!(photo_set, leaf, envir);

    return nothing
end
