###############################################################################
#
# Temperature dependency parameter set
#
###############################################################################
#= AbstractTDParameterSet type tree
---> ArrheniusTD
---> ArrheniusPeakTD
=#
"""
    type AbstractTDParameterSet{FT}

Hierarchy of the `AbstractTDParameterSet`:
- [`ArrheniusTD`](@ref)
- [`ArrheniusPeakTD`](@ref)
"""
abstract type AbstractTDParameterSet{FT} end




"""
    struct ArrheniusTD{FT}

An [`AbstractTDParameterSet`](@ref) type struct using
```math
corr = \\exp \\left( \\dfrac{ΔHa}{R T_0} - \\dfrac{ΔHa}{R T_1} \\right)
```

# Fields
$(DocStringExtensions.FIELDS)
"""
struct ArrheniusTD{FT<:AbstractFloat} <: AbstractTDParameterSet{FT}
    "Uncorrected value at 298.15 K"
    VAL_25     ::FT
    "Ratio between ΔHa and R `[K]`"
    ΔHa_to_R   ::FT
    "Ratio between ΔHa and R*K_25"
    ΔHa_to_RT25::FT
end




"""
    struct ArrheniusPeakTD{FT}

An [`AbstractTDParameterSet`](@ref) type struct using
```math
corr = \\exp \\left( \\dfrac{ΔHa}{R T_0} - \\dfrac{ΔHa}{R T_1} \\right)
       \\cdot
       \\dfrac{ 1 + \\exp \\left( \\dfrac{S_v T_0 - H_d}{R T_0} \\right) }
              { 1 + \\exp \\left( \\dfrac{S_v T_1 - H_d}{R T_1} \\right) }
```

# Fields
$(DocStringExtensions.FIELDS)
"""
struct ArrheniusPeakTD{FT<:AbstractFloat} <: AbstractTDParameterSet{FT}
    "Ratio between ΔHa and R*K_25"
    ΔHa_to_RT25::FT
    "Ratio between ΔHd and R"
    ΔHd_to_R   ::FT
    "Ratio between ΔSv and R"
    ΔSv_to_R   ::FT
    "Correction factor C = 1 + exp( Sv/R + Hd/(RT0) )"
    C          ::FT
end








###############################################################################
#
# Leaf fluorescence-related parameter set
#
###############################################################################
#= AbstractFluoModelParaSet type tree
---> FluoParaSet
=#
"""
    type AbstractFluoModelParaSet

Hierarchy of the `AbstractFluoModelParaSet`:
- [`FluoParaSet`](@ref)
"""
abstract type AbstractFluoModelParaSet{FT} end




"""
    struct FluoParaSet{FT}

A `AbstractFluoModelParaSet` type paramter set.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct FluoParaSet{FT<:AbstractFloat} <: AbstractFluoModelParaSet{FT}
    "Fluorescence model coefficient"
    Kn1::FT
    "Fluorescence model coefficient"
    Kn2::FT
    "Fluorescence model coefficient"
    Kn3::FT
end








###############################################################################
#
# Photosynthesis model parameter set -- temperature dependencies and etc
#
###############################################################################
#= AbstractPhotoModelParaSet type tree
---> C3ParaSet
---> C4ParaSet
=#
"""
    type AbstractPhotoModelParaSet

Hierarchy of the `AbstractPhotoModelParaSet`:
- [`C3ParaSet`](@ref)
- [`C4ParaSet`](@ref)
"""
abstract type AbstractPhotoModelParaSet{FT} end




"""
    struct C3Paraset{FT}

Parameter sets for C3 photosynthesis.

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct C3ParaSet{FT<:AbstractFloat} <: AbstractPhotoModelParaSet{FT}
    "Jmax temperature dependency"
    JT ::AbstractTDParameterSet{FT}
    "Kc temperature dependency"
    KcT::AbstractTDParameterSet{FT}
    "Ko temperature dependency"
    KoT::AbstractTDParameterSet{FT}
    "Respiration temperature dependency"
    ReT::AbstractTDParameterSet{FT}
    "Vcmax temperature dependency"
    VcT::AbstractTDParameterSet{FT}
    "Γ_star temperature dependency"
    ΓsT::AbstractTDParameterSet{FT}
    "Fluorescence model"
    Flu::AbstractFluoModelParaSet{FT}
    "Vcmax25 and respiration correlation"
    VR   ::FT
    "Coefficient 4.0/4.5 for NADPH/ATP requirement stochiometry, respectively"
    Eff_1::FT
    "Coefficient 8.0/10.5 for NADPH/ATP requirement stochiometry, respectively"
    Eff_2::FT
end




"""
    struct C4ParaSet{FT}

Parameter sets for C3 photosynthesis.

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct C4ParaSet{FT<:AbstractFloat} <: AbstractPhotoModelParaSet{FT}
    "Kpep temperature dependency"
    KpT::AbstractTDParameterSet{FT}
    "Respiration temperature dependency"
    ReT::AbstractTDParameterSet{FT}
    "Vcmax temperature dependency"
    VcT::AbstractTDParameterSet{FT}
    "Vpmax temperature dependency"
    VpT::AbstractTDParameterSet{FT}
    "Fluorescence model"
    Flu::AbstractFluoModelParaSet{FT}
    "Vcmax25 and respiration correlation"
    VR::FT
end








###############################################################################
#
# Leaf parameters container
#
###############################################################################
"""
    struct Leaf{FT}

Struct to store leaf information.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct Leaf{FT<:AbstractFloat}
    # Temperature related
    "Temperature `[K]`"
    T::FT = FT(298.15)

    # Photosynthesis system
    "Rate constant for thermal dissipation"
    Kd       ::FT = FT(0.85)
    "Rate constant for fluorescence (const)"
    Kf       ::FT = FT(0.05)
    "NPQ rate constant (initially zero)"
    Kn       ::FT = FT(0)
    "Rate constant for photochemistry (all reaction centers open)"
    Kp       ::FT = FT(4.0)
    "max PSII yield (Kn=0, all RC open)"
    maxPSII  ::FT = Kp/(Kp+Kf+Kd)
    "Fraction of absorbed light used by PSII ETR"
    PSII_frac::FT = FT(0.5)

    # CO₂ pressures
    "Leaf internal CO₂ partial pressure `[Pa]`"
    p_i  ::FT = FT(10)
    "Leaf surface CO₂ partial pressure `[Pa]`"
    p_s  ::FT = FT(40)
    "Saturation H₂O vapor pressure `[Pa]`"
    p_sat::FT = saturation_vapor_pressure(T)
    "Leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_bc ::FT = FT(3/1.35)
    "Leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_lc ::FT = FT(0.01)

    # Photosynthesis related
    "RubisCO limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Ac     ::FT = FT(0)
    "Light limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Aj     ::FT = FT(0)
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Ag     ::FT = FT(0)
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    An     ::FT = FT(0)
    "Product limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Ap     ::FT = FT(0)
    "Electron transport `[μmol m⁻² s⁻¹]`"
    J      ::FT = FT(0)
    "Potential Electron Transport Rate `[μmol m⁻² s⁻¹]`"
    J_pot  ::FT = FT(0)
    "Maximal electron transport rate `[μmol m⁻² s⁻¹]`"
    Jmax   ::FT = FT(120)
    "Maximal electron transport rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    Jmax25 ::FT = FT(120)
    "RubisCO coefficient Kc `[Pa]`"
    Kc     ::FT = FT(0)
    "RubisCO coefficient Ko `[Pa]`"
    Ko     ::FT = FT(0)
    "PEP coefficient Ko `[Pa]`"
    Kpep   ::FT = FT(0)
    "Michaelis-Menten's coefficient `[Pa]`"
    Km     ::FT = FT(0)
    "Respiration rate `[μmol m⁻² s⁻¹]`"
    Rd     ::FT = FT(1)
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    Rd25   ::FT = FT(1)
    "Maximal carboxylation rate `[μmol m⁻² s⁻¹]`"
    Vcmax  ::FT = FT(60)
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    Vcmax25::FT = FT(60)
    "Maximal PEP carboxylation rate `[μmol m⁻² s⁻¹]`"
    Vpmax  ::FT = FT(120)
    "Maximal PEP carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    Vpmax25::FT = FT(120)
    "CO₂ compensation point with the absence of Rd `[Pa]`"
    Γ_star ::FT = FT(0)

    # Fluorescence related
    "Total efficiency, incl. photorespiration `[mol CO₂ mol⁻¹ e-]`"
    CO₂_per_electron::FT = FT(1/6)
    "dark adapted yield (`Kp=0`)"
    Fm              ::FT = FT(0)
    "light adapted yield (`Kp=0`)"
    Fm′             ::FT = FT(0)
    "dark-adapted fluorescence yield (`Kp=max`)"
    Fo              ::FT = FT(0)
    "light-adapted fluorescence yield in the dark (`Kp=max`)"
    Fo′             ::FT = FT(0)
    "Actual electron transport rate `[μmol m⁻² s⁻¹]`"
    Ja              ::FT = FT(0)
    "Non-Photochemical quenching "
    NPQ             ::FT = FT(0)
    "Photochemical quenching"
    qQ              ::FT = FT(0)
    "energy quenching"
    qE              ::FT = FT(0)
    "PSII yield"
    φ               ::FT = FT(0)
    "Steady-state (light-adapted) yield (aka Fs)"
    ϕs              ::FT = FT(0)

    # Environment related
    "Absorbed photosynthetic active radiation `[μmol m⁻² s⁻¹]`"
    APAR::FT = FT(100)
end








###############################################################################
#
# Environmental conditions container
#
###############################################################################
"""
    struct AirLayer{FT}

Struct to store environmental conditions in each air layer corresponds to one
    canopy layer.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct AirLayer{FT<:AbstractFloat}
    "Air temperature `[K]`"
    t_air::FT = 298.15

    # pressures
    "Atmospheric CO₂ partial pressure `[Pa]`"
    p_a  ::FT = FT(41.0)
    "Atmospheric pressure `[Pa]`"
    p_atm::FT = FT(101325.0)
    "Atmospheric vapor pressure `[Pa]`"
    p_H₂O::FT = FT(1500.0)
    "Atmospheric O₂ partial pressure `[Pa]`"
    p_O₂ ::FT = FT(101325.0 * 0.209)
    "Saturation vapor pressure `[Pa]`"
    p_sat::FT = saturation_vapor_pressure(t_air)
    "Relatiev humidity"
    RH   ::FT = p_H₂O / p_sat
end
