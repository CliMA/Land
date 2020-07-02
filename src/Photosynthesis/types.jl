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
    type AbstractTDParameterSet

Hierarchy of the `AbstractTDParameterSet`:
- [`ArrheniusTD`](@ref)
- [`ArrheniusPeakTD`](@ref)
"""
abstract type AbstractTDParameterSet end




"""
    struct ArrheniusTD{FT}

An [`AbstractTDParameterSet`](@ref) type struct using
```math
corr = \\exp \\left( \\dfrac{ΔHa}{R T_0} - \\dfrac{ΔHa}{R T_1} \\right)
```

# Fields
$(DocStringExtensions.FIELDS)
"""
struct ArrheniusTD{FT<:AbstractFloat} <: AbstractTDParameterSet
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
struct ArrheniusPeakTD{FT<:AbstractFloat} <: AbstractTDParameterSet
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
abstract type AbstractStomatalModel end




"""
    type EmpiricalStomatalModel

Hierarchy of the `EmpiricalStomatalModel`:
- [`ESMBallBerry`](@ref)
- [`ESMGentine`](@ref)
- [`ESMLeuning`](@ref)
- [`ESMMedlyn`](@ref)
"""
abstract type EmpiricalStomatalModel <: AbstractStomatalModel end




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
Base.@kwdef mutable struct ESMBallBerry{FT<:AbstractFloat} <: EmpiricalStomatalModel
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
Base.@kwdef mutable struct ESMGentine{FT<:AbstractFloat} <: EmpiricalStomatalModel
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
Base.@kwdef mutable struct ESMLeuning{FT<:AbstractFloat} <: EmpiricalStomatalModel
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
Base.@kwdef mutable struct ESMMedlyn{FT<:AbstractFloat} <: EmpiricalStomatalModel
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
abstract type OptimizationStomatalModel <: AbstractStomatalModel end




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
struct OSMEller <: OptimizationStomatalModel end




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
struct OSMSperry <: OptimizationStomatalModel end




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
struct OSMWang <: OptimizationStomatalModel end




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
Base.@kwdef mutable struct OSMWAP{FT<:AbstractFloat} <: OptimizationStomatalModel
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
Base.@kwdef mutable struct OSMWAPMod{FT<:AbstractFloat} <: OptimizationStomatalModel
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
abstract type AbstractFluoModelParaSet end




"""
    struct FluoParaSet{FT}

A `AbstractFluoModelParaSet` type paramter set.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct FluoParaSet{FT<:AbstractFloat} <: AbstractFluoModelParaSet
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
abstract type AbstractPhotoModelParaSet end




"""
    struct C3Paraset{FT}

Parameter sets for C3 photosynthesis.

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct C3ParaSet{FT<:AbstractFloat} <: AbstractPhotoModelParaSet
    "Jmax temperature dependency"
    JT ::AbstractTDParameterSet
    "Kc temperature dependency"
    KcT::AbstractTDParameterSet
    "Ko temperature dependency"
    KoT::AbstractTDParameterSet
    "Respiration temperature dependency"
    ReT::AbstractTDParameterSet
    "Vcmax temperature dependency"
    VcT::AbstractTDParameterSet
    "Γ_star temperature dependency"
    ΓsT::AbstractTDParameterSet
    "Colimitation mode"
    Col::AbstractColimitation
    "Fluorescence model"
    Flu::AbstractFluoModelParaSet
    "Stomatal model scheme"
    Sto::AbstractStomatalModel
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
mutable struct C4ParaSet{FT<:AbstractFloat} <: AbstractPhotoModelParaSet
    "Kpep temperature dependency"
    KpT::AbstractTDParameterSet
    "Respiration temperature dependency"
    ReT::AbstractTDParameterSet
    "Vcmax temperature dependency"
    VcT::AbstractTDParameterSet
    "Vpmax temperature dependency"
    VpT::AbstractTDParameterSet
    "Colimitation mode"
    Col::AbstractColimitation
    "Fluorescence model"
    Flu::AbstractFluoModelParaSet
    "Stomatal model scheme"
    Sto::AbstractStomatalModel
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
    "Sensible Heat Flux `[W m⁻²]`"
    H    ::FT = FT(0)
    "Latent Heat Flux `[W m⁻²]`"
    LE   ::FT = FT(0)
    "Latent Heat of evaporation `[J mol⁻¹]`"
    LV   ::FT = latent_heat_vapor(K_25) * 1000 / MOLMASS_WATER
    "Net Radiation Balance `[W m⁻²]`"
    Rn   ::FT = FT(0)
    "Temperature `[K]`"
    T    ::FT = FT(K_25)
    "Old temperature `[K]`"
    T_old::FT = FT(0)
    "Leaf width `[m]`"
    width::FT = FT(0.05)

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

    # Diffusive conductances
    "Boundary layer conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_bc   ::FT = FT(3/1.35)
    "Boundary layer conductance to heat `[mol m⁻² s⁻¹]`"
    g_bh   ::FT = FT(3)
    "Boundary layer conductance to H₂O `[mol m⁻² s⁻¹]`"
    g_bw   ::FT = FT(3)
    "Gias correction constant"
    g_ias_c::FT = FT(0)
    "Gias correction exponent"
    g_ias_e::FT = FT(0.3)
    "Leaf diffusive conductance to water CO₂ `[mol m⁻² s⁻¹]`"
    g_lc   ::FT = FT(0.08032)
    "Leaf diffusive conductance to water H₂O `[mol m⁻² s⁻¹]`"
    g_lw   ::FT = FT(0.1519)
    "Mesophyll conductance for CO₂ `[mol m⁻² s⁻¹]`"
    g_m    ::FT = FT(0.5)
    "Maximal leaf diffusive conductance `[mol m⁻² s⁻¹]`"
    g_max  ::FT = FT(0.8)
    "Maximal leaf diffusive conductance at 298.15 K `[mol m⁻² s⁻¹]`"
    g_max25::FT = FT(0.8)
    "Minimal leaf diffusive conductance `[mol m⁻² s⁻¹]`"
    g_min  ::FT = FT(0.025)
    "Minimal leaf diffusive conductance at 298.15 K `[mol m⁻² s⁻¹]`"
    g_min25::FT = FT(0.025)
    "Stomatal conductance to water CO₂ `[mol m⁻² s⁻¹]`"
    g_sc   ::FT = FT(0.1)
    "Stomatal conductance to water H₂O `[mol m⁻² s⁻¹]`"
    g_sw   ::FT = FT(0.16)

    # CO₂ pressures
    "Leaf internal CO₂ partial pressure `[Pa]`"
    p_i  ::FT = FT(10)
    "Leaf surface CO₂ partial pressure `[Pa]`"
    p_s  ::FT = FT(40)
    "Leaf saturation vapor pressure `[Pa]`"
    p_sat::FT = saturation_vapor_pressure(T)

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

    # Leaf hydraulics system and stomtal optimization related
    "Leaf hydraulic system"
    hs::LeafHydraulics = LeafHydraulics{FT}()
    "Maximal photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_max ::FT = FT(0)
    "Flow rate `[mol m⁻² s⁻¹]`"
    e     ::FT = FT(0)
    "Critical flow rate `[mol m⁻² s⁻¹]`"
    ec    ::FT = FT(2e-9)
    "Maximal hydraulic conductance ratio"
    kr_max::FT = FT(1)
    "Base xylem pressre `[MPa]`"
    p_ups ::FT = FT(0)
    "Base xylem pressre memory `[MPa]`"
    p_old ::FT = FT(1)
end




"""
    struct Leaves{FT}

Struct to store leaf information.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct Leaves{FT<:AbstractFloat}
    # Number of Leaves per canopy layer
    n_leaf::Int = 325

    # Temperature related, different for each leaf
    "Sensible Heat Flux `[W m⁻²]`"
    H ::Array{FT,1} = zeros(FT, n_leaf)
    "Latent Heat Flux `[W m⁻²]`"
    LE::Array{FT,1} = zeros(FT, n_leaf)
    "Net Radiation Balance `[W m⁻²]`"
    Rn::Array{FT,1} = zeros(FT, n_leaf)

    # Tempearture related, same for all leaves
    "Latent Heat of evaporation `[J mol⁻¹]`"
    LV   ::FT = latent_heat_vapor(K_25) * 1000 / MOLMASS_WATER
    "Temperature `[K]`"
    T    ::FT = FT(K_25)
    "Old temperature `[K]`"
    T_old::FT = FT(0)
    "Leaf width `[m]`"
    width::FT = FT(0.05)

    # Photosynthesis system, different for each leaf
    "NPQ rate constant (initially zero)"
    Kn     ::Array{FT,1} = zeros(FT, n_leaf)
    "Rate constant for photochemistry (all reaction centers open)"
    Kp     ::Array{FT,1} = zeros(FT, n_leaf) .+ 4

    # Photosynthesis system, same for all leaves
    "Rate constant for thermal dissipation"
    Kd       ::FT = FT(0.85)
    "Rate constant for fluorescence (const)"
    Kf       ::FT = FT(0.05)
    "Maximal rate constant for photochemistry (all reaction centers open)"
    Kpmax    ::FT = FT(4)
    "max PSII yield (Kn=0, all RC open)"
    maxPSII  ::FT = Kpmax / (Kpmax + Kf +Kd)
    "Fraction of absorbed light used by PSII ETR"
    PSII_frac::FT = FT(0.5)

    # Diffusive conductances, same for all leaves
    "Boundary layer conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_bc::Array{FT,1} = zeros(FT, n_leaf) .+ FT(3/1.35)
    "Boundary layer conductance to heat `[mol m⁻² s⁻¹]`"
    g_bh::Array{FT,1} = zeros(FT, n_leaf) .+ FT(3)
    "Boundary layer conductance to H₂O `[mol m⁻² s⁻¹]`"
    g_bw::Array{FT,1} = zeros(FT, n_leaf) .+ FT(3)
    "Leaf diffusive conductance to water CO₂ `[mol m⁻² s⁻¹]`"
    g_lc::Array{FT,1} = zeros(FT, n_leaf) .+ FT(0.08032)
    "Leaf diffusive conductance to water H₂O `[mol m⁻² s⁻¹]`"
    g_lw::Array{FT,1} = zeros(FT, n_leaf) .+ FT(0.1519)
    "Mesophyll conductance for CO₂ `[mol m⁻² s⁻¹]`"
    g_m ::Array{FT,1} = zeros(FT, n_leaf) .+ FT(0.5)
    "Stomatal conductance to water CO₂ `[mol m⁻² s⁻¹]`"
    g_sc::Array{FT,1} = zeros(FT, n_leaf) .+ FT(0.1)
    "Stomatal conductance to water H₂O `[mol m⁻² s⁻¹]`"
    g_sw::Array{FT,1} = zeros(FT, n_leaf) .+ FT(0.16)

    # Diffusive conductances, same for all leaves
    "Gias correction constant"
    g_ias_c::FT = FT(0)
    "Gias correction exponent"
    g_ias_e::FT = FT(0.3)
    "Maximal leaf diffusive conductance `[mol m⁻² s⁻¹]`"
    g_max  ::FT = FT(0.8)
    "Maximal leaf diffusive conductance at 298.15 K `[mol m⁻² s⁻¹]`"
    g_max25::FT = FT(0.8)
    "Minimal leaf diffusive conductance `[mol m⁻² s⁻¹]`"
    g_min  ::FT = FT(0.025)
    "Minimal leaf diffusive conductance at 298.15 K `[mol m⁻² s⁻¹]`"
    g_min25::FT = FT(0.025)

    # CO₂ and H₂O pressures, different for each leaf
    "Leaf internal CO₂ partial pressure `[Pa]`"
    p_i::Array{FT,1} = zeros(FT, n_leaf) .+ FT(10)
    "Leaf surface CO₂ partial pressure `[Pa]`"
    p_s::Array{FT,1} = zeros(FT, n_leaf) .+ FT(40)

    # CO₂ and H₂O pressures, same for all leaves
    "Leaf saturation vapor pressure `[Pa]`"
    p_sat::FT = saturation_vapor_pressure(T)

    # Photosynthesis related, different for each leaf
    "RubisCO limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Ac   ::Array{FT,1} = zeros(FT, n_leaf)
    "Light limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Aj   ::Array{FT,1} = zeros(FT, n_leaf)
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Ag   ::Array{FT,1} = zeros(FT, n_leaf)
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    An   ::Array{FT,1} = zeros(FT, n_leaf)
    "Product limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Ap   ::Array{FT,1} = zeros(FT, n_leaf)
    "Electron transport `[μmol m⁻² s⁻¹]`"
    J    ::Array{FT,1} = zeros(FT, n_leaf)
    "Potential Electron Transport Rate `[μmol m⁻² s⁻¹]`"
    J_pot::Array{FT,1} = zeros(FT, n_leaf)

    # Photosynthesis related, same for all leaves
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

    # Fluorescence related, different for each leaf
    "Total efficiency, incl. photorespiration `[mol CO₂ mol⁻¹ e-]`"
    CO₂_per_electron::Array{FT,1} = zeros(FT, n_leaf) .+ FT(1/6)
    "light adapted yield (`Kp=0`)"
    Fm′             ::Array{FT,1} = zeros(FT, n_leaf)
    "light-adapted fluorescence yield in the dark (`Kp=max`)"
    Fo′             ::Array{FT,1} = zeros(FT, n_leaf)
    "Actual electron transport rate `[μmol m⁻² s⁻¹]`"
    Ja              ::Array{FT,1} = zeros(FT, n_leaf)
    "Non-Photochemical quenching "
    NPQ             ::Array{FT,1} = zeros(FT, n_leaf)
    "Photochemical quenching"
    qQ              ::Array{FT,1} = zeros(FT, n_leaf)
    "energy quenching"
    qE              ::Array{FT,1} = zeros(FT, n_leaf)
    "PSII yield"
    φ               ::Array{FT,1} = zeros(FT, n_leaf)
    "Steady-state (light-adapted) yield (aka Fs)"
    ϕs              ::Array{FT,1} = zeros(FT, n_leaf)

    # Fluorescence related, same for all leaves
    "dark adapted yield (`Kp=0`)"
    Fm::FT = FT(0)
    "dark-adapted fluorescence yield (`Kp=max`)"
    Fo::FT = FT(0)

    # Environment related, different for each leaf
    "Absorbed photosynthetic active radiation `[μmol m⁻² s⁻¹]`"
    APAR::Array{FT,1} = zeros(FT, n_leaf) .+ 100

    # Leaf hydraulics system
    "Leaf hydraulic system"
    hs::LeafHydraulics = LeafHydraulics{FT}()

    # Stomtal optimization related, different for each leaf
    "Maximal photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_max::Array{FT,1} = zeros(FT, n_leaf)

    # Stomtal optimization related, same for all leaves
    "Flow rate `[mol m⁻² s⁻¹]`"
    e     ::FT = FT(0)
    "Critical flow rate `[mol m⁻² s⁻¹]`"
    ec    ::FT = FT(2e-9)
    "Maximal hydraulic conductance ratio"
    kr_max::FT = FT(1)
    "Base xylem pressre `[MPa]`"
    p_ups ::FT = FT(0)
    "Base xylem pressre memory `[MPa]`"
    p_old ::FT = FT(1)
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
    t_air::FT = K_25

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
