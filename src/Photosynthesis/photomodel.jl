###############################################################################
#
# Arrhenius corrections
#
###############################################################################
"""
    arrhenius_correction(td_set::ArrheniusTD{FT}    , T::FT)
    arrhenius_correction(td_set::ArrheniusPeakTD{FT}, T::FT)

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
            td_set::AbstractTDParameterSet,
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
            td_set::AbstractTDParameterSet,
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
            td_set::AbstractTDParameterSet,
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
            td_set::AbstractTDParameterSet,
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
            td_set::AbstractTDParameterSet,
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
    @unpack APAR, maxPSII, Jmax, θ_j, PSII_frac = leaf

    _Jp = PSII_frac * maxPSII * APAR
    _b  = _Jp + Jmax
    _c  = _Jp * Jmax
    _J  = ( _b - sqrt(_b^2 - 4*θ_j*_c) ) / (2*θ_j)

    leaf.J_pot = _Jp
    leaf.J     = _J

    return nothing
end

function leaf_ETR!(
            photo_set::C4ParaSet{FT},
            leaf::Leaf{FT}
            ) where {FT<:AbstractFloat}
    @unpack APAR, maxPSII, PSII_frac = leaf

    _Jp = PSII_frac * maxPSII * APAR

    leaf.J_pot = _Jp
    leaf.J     = _Jp

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
    leaf.Aj               = leaf.J * leaf.CO₂_per_electron;

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
# Compute the colimited photosynthetic rate
#
###############################################################################
"""
    photosynthesis_colimit(colim::MinColimit, A...)
    photosynthesis_colimit(colim::CurvedColimit{FT}, A...)

Return the colimited photosynthetic rate, given
-`colim` [`AbstractColimitation`](@ref) type co-limitation
-`A` Tuple of photosynthetic rates

Multiple options were added to save time
"""
function photosynthesis_colimit(
            colim::MinColimit,
            A...)
    return min(A...)
end

function photosynthesis_colimit(
            colim::CurvedColimit{FT},
            A1::FT,
            A2::FT
            ) where {FT<:AbstractFloat}
    a = lower_quadratic(colim.curvature, -(A1 + A2), A2 * A2);

    return isnan(a) ? min(A1,A2) : a
end

function photosynthesis_colimit(
            colim::CurvedColimit{FT},
            A1::FT,
            A2::FT,
            A3::FT
            ) where {FT<:AbstractFloat}
    a = lower_quadratic(colim.curvature, -(A1 + A2), A1 * A2);
    a = lower_quadratic(colim.curvature, -(a  + A3), a  * A3);

    return isnan(a) ? min(A1,A2,A3) : a
end

function photosynthesis_colimit(
            colim::CurvedColimit{FT},
            A...
            ) where {FT<:AbstractFloat}
    a = A[1];
    for j=2:length(A)
        a = lower_quadratic(colim.curvature, -(a + A[j]), a * A[j]);
    end

    return isnan(a) ? minimum(A) : a
end








###############################################################################
#
# Update the fluorescence in the leaf
#
###############################################################################
"""
    leaf_fluorescence!(fluo_set::FluoParaSet{FT}, leaf::Leaf{FT})

Compute fluorescence yield, Kn, and Kp for leaf, given
- `fluo_set` [`FluoParaSet`](@ref) type parameter set
- `leaf` [`Leaf`](@ref) struct
"""
function leaf_fluorescence!(
            fluo_set::FluoParaSet{FT},
            leaf::Leaf{FT}
            ) where {FT<:AbstractFloat}
    @unpack Ag, Kd, Kf, maxPSII = leaf;
    @unpack Kn1, Kn2, Kn3 = fluo_set;

    # Actual effective ETR:
    leaf.Ja  = max(0, Ag / leaf.CO₂_per_electron);

    # Effective photochemical yield:
    if leaf.Ja <= 0
        _φ   = maxPSII;
    else
        _φ   = maxPSII*leaf.Ja/leaf.J_pot;
    end

    leaf.φ   = min(1/maxPSII, _φ);
    # degree of light saturation: 'x' (van der Tol e.Ap. 2014)
    x        = max(0,  1-leaf.φ/maxPSII);
    
    # Max PSII rate constant
    Kp_max   = FT(4.0);

    x_alpha  = exp(log(x)*Kn2);
    #println(x_alpha)
    
    leaf.Kn  = Kn1 * (1+Kn3)* x_alpha/(Kn3 + x_alpha);
    leaf.Kp  = max(0,-leaf.φ*(Kf+Kd+leaf.Kn)/(leaf.φ-1));

    leaf.Fo  = Kf/(Kf+Kp_max+Kd        );
    leaf.Fo′ = Kf/(Kf+Kp_max+Kd+leaf.Kn);
    leaf.Fm  = Kf/(Kf       +Kd        );
    leaf.Fm′ = Kf/(Kf       +Kd+leaf.Kn);
    leaf.ϕs  = leaf.Fm′*(1-leaf.φ);
    # leaf.eta  = leaf.ϕs/leaf.Fo
    # don't need this anymore
    # better to use ϕs directly for SIF as Fo is not always fqe=0.01
    leaf.qQ  = 1-(leaf.ϕs-leaf.Fo′)/(leaf.Fm-leaf.Fo′);
    leaf.qE  = 1-(leaf.Fm-leaf.Fo′)/(leaf.Fm′-leaf.Fo);
    leaf.NPQ = leaf.Kn/(Kf+Kd);

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
            leaf::Leaf{FT},
            envir::AirLayer{FT})

Update the temperature dependent photosynthesis only, given
- `photo_set` [`AbstractPhotoModelParaSet`](@ref) type parameter set
- `leaf` [`Leaf`](@ref) type struct
- `envir` [`AirLayer`](@ref) type struct
"""
function leaf_temperature_dependence!(
            photo_set::C3ParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT}
            ) where {FT<:AbstractFloat}
    leaf.g_max      = leaf.g_max25 * relative_diffusive_coefficient(leaf.T);
    leaf.g_min      = leaf.g_min25 * relative_diffusive_coefficient(leaf.T);
    leaf.LV         = latent_heat_vapor(leaf.T) * 1000 / FT(MOLMASS_WATER);
    leaf.p_sat      = saturation_vapor_pressure(leaf.T);
    (leaf.hs).f_st  = relative_surface_tension(leaf.T);
    (leaf.hs).f_vis = relative_viscosity(leaf.T);

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
    leaf.g_max      = leaf.g_max25 * relative_diffusive_coefficient(leaf.T);
    leaf.g_min      = leaf.g_min25 * relative_diffusive_coefficient(leaf.T);
    leaf.LV         = latent_heat_vapor(leaf.T) * 1000 / FT(MOLMASS_WATER);
    leaf.p_sat      = saturation_vapor_pressure(leaf.T);
    (leaf.hs).f_st  = relative_surface_tension(leaf.T);
    (leaf.hs).f_vis = relative_viscosity(leaf.T);

    leaf_rd!(photo_set.ReT, leaf);
    leaf_vcmax!(photo_set.VcT, leaf);
    leaf_vpmax!(photo_set.VpT, leaf);
    leaf_kpep!(photo_set.KpT, leaf);

    return nothing
end








###############################################################################
#
# Calculate photosynthesis from CO₂ partial pressure
#
###############################################################################
"""
    leaf_photo_from_pi!(
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT})

Compute leaf photosynthetic rates, given
- `photo_set` [`AbstractPhotoModelParaSet`](@ref) type parameter set
- `leaf` [`Leaf`](@ref) type struct
- `envir` [`AirLayer`](@ref) type struct

The C3 photosynthesis model is from Farquhar et al. (1980) "A biochemical model
    of photosynthetic CO₂ assimilation in leaves of C3 species."

The C4 photosynthesis model is adapted from Collatz et al. (1992) "Coupled
    photosynthesis-stomatal conductance model for leaves of C4 plants."
"""
function leaf_photo_from_pi!(
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT}
            ) where {FT<:AbstractFloat}
    # update TD only if T changes
    if leaf.T != leaf.T_old
        leaf_temperature_dependence!(photo_set, leaf, envir);
        leaf.T_old = leaf.T;
    end

    leaf_ETR!(photo_set, leaf);
    light_limited_rate!(photo_set, leaf);
    rubisco_limited_rate!(photo_set, leaf);
    product_limited_rate!(photo_set, leaf);
    leaf.Ag = photosynthesis_colimit(photo_set.Col, leaf.Ac, leaf.Aj, leaf.Ap);
    leaf.An = leaf.Ag - leaf.Rd;

    return nothing
end








###############################################################################
#
# Calculate photosynthesis from leaf diffusive conductance for CO₂
#
###############################################################################
"""
    leaf_photo_from_glc!(
            photo_set::AbstractPhotoModelParaSet,
            leaf::Leaf{FT},
            envir::AirLayer{FT})

Update leaf photosynthetic rates from a known leaf diffusive conductance, given
- `photo_set` [`AbstractPhotoModelParaSet`](@ref) type parameter set
- `leaf` [`Leaf`](@ref) type struct
- `envir` [`AirLayer`](@ref) type struct
"""
function leaf_photo_from_glc!(
            photo_set::C3ParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT}
            ) where {FT<:AbstractFloat}
    # update TD only if T changes
    if leaf.T != leaf.T_old
        leaf_temperature_dependence!(photo_set, leaf, envir);
        leaf.T_old = leaf.T;
    end

    leaf_ETR!(photo_set, leaf);
    light_limited_rate_glc!(photo_set, leaf, envir);
    rubisco_limited_rate_glc!(photo_set, leaf, envir);
    product_limited_rate!(photo_set, leaf);
    leaf.Ag  = photosynthesis_colimit(photo_set.Col, leaf.Ac, leaf.Aj, leaf.Ap);
    leaf.An  = leaf.Ag - leaf.Rd;
    leaf.p_i = envir.p_a - leaf.An*FT(1e-6) * envir.p_atm / leaf.g_lc;
    leaf.p_s = envir.p_a - leaf.An*FT(1e-6) * envir.p_atm / leaf.g_bc;

    return nothing
end

function leaf_photo_from_glc!(
            photo_set::C4ParaSet{FT},
            leaf::Leaf{FT},
            envir::AirLayer{FT}
            ) where {FT<:AbstractFloat}
    # update TD only if T changes
    if leaf.T != leaf.T_old
        leaf_temperature_dependence!(photo_set, leaf, envir);
        leaf.T_old = leaf.T;
    end

    leaf_ETR!(photo_set, leaf);
    light_limited_rate!(photo_set, leaf);
    rubisco_limited_rate!(photo_set, leaf);
    product_limited_rate_glc!(photo_set, leaf, envir);
    leaf.Ag  = photosynthesis_colimit(photo_set.Col, leaf.Ac, leaf.Aj, leaf.Ap);
    leaf.An  = leaf.Ag - leaf.Rd;
    leaf.p_i = envir.p_a - leaf.An*FT(1e-6) * envir.p_atm / leaf.g_lc;
    leaf.p_s = envir.p_a - leaf.An*FT(1e-6) * envir.p_atm / leaf.g_bc;

    return nothing
end
