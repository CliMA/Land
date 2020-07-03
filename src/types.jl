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
# Stomatal model scheme types
#
###############################################################################
#= EmpiricalStomatalModel type tree
AbstractStomatalModel
---> EmpiricalStomatalModel
    ---> ESMBallBerry    # use RH
    ---> ESMGentine      # use k_leaf
    ---> ESMLeuning      # Use VPD
    ---> ESMMedlyn       # use VPD
---> OptimizationStomatalModel
    ---> OSMEller        # Eller 2018 Model
    ---> OSMSperry       # Sperry 2017 model
    ---> OSMWang         # Wang 2020 model
    ---> OSMWAP          # Wolf-Anderegg-Pacala model
    ---> OSMWAP          # Modified Wolf-Anderegg-Pacala model
=#
"""
    type AbstractStomatalModel

Hierarchy of the `AbstractStomatalModel`:
- [`EmpiricalStomatalModel`](@ref)
- [`OptimizationStomatalModel`](@ref)
"""
abstract type AbstractStomatalModel{FT} end




"""
    type EmpiricalStomatalModel

Hierarchy of the `EmpiricalStomatalModel`:
- [`ESMBallBerry`](@ref)
- [`ESMGentine`](@ref)
- [`ESMLeuning`](@ref)
- [`ESMMedlyn`](@ref)
"""
abstract type EmpiricalStomatalModel{FT} <: AbstractStomatalModel{FT} end




"""
    struct ESMBallBerry{FT}

An empirical model parameter set type for Ball-Berry type model.
The equation used for Ball-Berry type model is
```math
gs = g0 + g1 ⋅ RH ⋅ \\dfrac{A}{Cs}
```

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct ESMBallBerry{FT<:AbstractFloat} <: EmpiricalStomatalModel{FT}
    "minimal stomatal conductance g0 `[mol m⁻² s⁻¹]`"
    g0::FT = FT(0.025)
    "slope of conductance-photosynthesis correlation `[unitless]`"
    g1::FT = FT(9.0  )
end




"""
    struct ESMGentine{FT}

An empirical model parameter set type for Gentine type model.
The equation used for Gentine type model is
```math
gs = g0 + g1 ⋅ \\dfrac{k_{leaf}}{k_{max}} ⋅ \\dfrac{A}{Ca}.
```
Note it that the Gentine model does not require for a `β` function to tune the
    soil drought response, but the use of `k_leaf` also does not permit
    post-drought stomatal response unless `k_leaf` can be recovered.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct ESMGentine{FT<:AbstractFloat} <: EmpiricalStomatalModel{FT}
    "minimal stomatal conductance g0 `[mol m⁻² s⁻¹]`"
    g0::FT = FT(0.025)
    "slope of conductance-photosynthesis correlation `[unitless]`"
    g1::FT = FT(9.0  )
end




"""
    struct ESMLeuning{FT}

An empirical model parameter set type for Leuning type model.
The equation used for Leuning type model is
```math
gs = g0 + g1 ⋅ \\dfrac{A}{Cs - Γ^{*}} ⋅ \\dfrac{1}{1 + \\dfrac{VPD}{d0}}
```

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct ESMLeuning{FT<:AbstractFloat} <: EmpiricalStomatalModel{FT}
    "minimal stomatal conductance g0 `[mol m⁻² s⁻¹]`"
    g0::FT = FT(0.025 )
    "slope of conductance-photosynthesis correlation `[unitless]`"
    g1::FT = FT(8.0   )
    "fitting parameter of d/d0 below the fraction, same unit as vpd `[Pa]`"
    d0::FT = FT(3000.0)
end




"""
    struct ESMMedlyn{FT}

An empirical model parameter set type for Medlyn type model.
The equation used in Medlyn type model is
```
gs = g0 + (1 + \\dfrac{g1}{\\sqrt{VPD}} ⋅ \\dfrac{A}{Ca}
```

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct ESMMedlyn{FT<:AbstractFloat} <: EmpiricalStomatalModel{FT}
    "minimal stomatal conductance g0 `[mol m⁻² s⁻¹]`"
    g0::FT = FT(0.025)
    "slope of conductance-photosynthesis correlation `[Pa⁽⁵⁾]`"
    g1::FT = FT(125.0)
end




"""
    type OptimizationStomatalModel

Hierarchy of the `OptimizationStomatalModel`:
- [`OSMEller`](@ref)
- [`OSMSperry`](@ref)
- [`OSMWang`](@ref)
- [`OSMWAP`](@ref)
- [`OSMWAPMod`](@ref)
"""
abstract type OptimizationStomatalModel{FT} <: AbstractStomatalModel{FT} end




"""
    struct OSMEller

An optimization model parameter set type for Eller model.
The equation used for Eller model is
```math
\\dfrac{∂Θ}{∂E} = -\\dfrac{∂K}{∂E} ⋅ \\dfrac{A}{K}
```
where K is ``\\dfrac{∂E}{∂P}``.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct OSMEller{FT} <: OptimizationStomatalModel{FT} end




"""
    struct OSMSperry

An optimization model parameter set type for Sperry model.
The equation used for Sperry model is
```math
\\dfrac{∂Θ}{∂E} = -\\dfrac{∂K}{∂E} ⋅ \\dfrac{A_{max}}{K_{max}}
```
where K is ``\\dfrac{∂E}{∂P}``.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct OSMSperry{FT} <: OptimizationStomatalModel{FT} end




"""
    struct OSMWang

An optimization model parameter set type for Eller type model.
The equation used for Wang model is
```math
\\dfrac{∂Θ}{∂E} = \\dfrac{A}{E_{crit} - E}
```

# Fields
$(DocStringExtensions.FIELDS)
"""
struct OSMWang{FT} <: OptimizationStomatalModel{FT} end




"""
    struct OSMWAP{FT}

An optimization model parameter set type for Wolf-Anderegg-Pacala type model.
The equation used for Wolf-Anderegg-Pacala model is
```math
\\dfrac{∂Θ}{∂E} = \\dfrac{2aP + b}{K}
```
where K is ∂P/∂E.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct OSMWAP{FT<:AbstractFloat} <: OptimizationStomatalModel{FT}
    "Quadratic equation parameter `[μmol m⁻² s⁻¹ MPa⁻²]`"
    a::FT = FT(0.5)
    "Quadratic equation parameter `[μmol m⁻² s⁻¹ MPa⁻¹]`"
    b::FT = FT(2.0)
end




"""
    struct OSMWAPMod{FT}

An optimization model parameter set type for Wolf-Anderegg-Pacala type model,
    modified by adding a photosynthesis component while set b and c = 0.
The equation used for modified Wolf-Anderegg-Pacala model is
```math
\\dfrac{∂Θ}{∂E} = \\dfrac{aAP}{K}
```
where P is absolute value of leaf xylem pressure.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct OSMWAPMod{FT<:AbstractFloat} <: OptimizationStomatalModel{FT}
    "Quadratic equation parameter `[mol mol⁻¹ MPa⁻¹]`"
    a::FT = FT(0.1)
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
    "Stomatal model scheme"
    Sto::AbstractStomatalModel{FT}
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
    "Stomatal model scheme"
    Sto::AbstractStomatalModel{FT}
    "Vcmax25 and respiration correlation"
    VR::FT
end








###############################################################################
#
# Leaf parameters container
#
###############################################################################
"""
    type AbstractLeaf

Hierarchy of AbstractLeaf
- [`Leaf`](@ref)
- [`Leaves`](@ref)
"""
abstract type AbstractLeaf{FT} end




"""
    struct Leaf{FT}

Struct to store leaf information.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct Leaf{FT<:AbstractFloat} <: AbstractLeaf{FT}
    # Temperature related
    "Temperature `[K]`"
    T    ::FT = FT(298.15)

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
    p_sat::FT = FT(3166.0)
    "Relatiev humidity"
    RH   ::FT = p_H₂O / p_sat
end
