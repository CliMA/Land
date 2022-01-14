###############################################################################
#
# Environmental conditions cache
#
###############################################################################
"""
    mutable struct AirLayer{FT}

Struct to store environmental conditions in each air layer corresponds to one
    canopy layer.

# Fields
$(TYPEDFIELDS)

---
# Examples
```julia
# create a layer of air surrouding the leaf
envir = AirLayer{Float32}();
```
"""
Base.@kwdef mutable struct AirLayer{FT<:AbstractFloat}
    # pressures
    "Atmospheric CO₂ partial pressure `[Pa]`"
    p_a::FT = 41
    "Atmospheric pressure `[Pa]`"
    p_atm::FT = 101325
    "Atmospheric O₂ partial pressure `[Pa]`"
    p_O₂::FT = p_atm * 0.209
end








###############################################################################
#
# Leaf parameters container
#
###############################################################################
"""
    mutable struct Leaf{FT}

Struct to store leaf information.

# Fields
$(TYPEDFIELDS)

---
# Examples
```julia
# create a leaf (both C3 and C4 competible)
leaf = Leaf{Float32}();
```
"""
Base.@kwdef mutable struct Leaf{FT<:AbstractFloat}
    # Temperature related
    "Temperature `[K]`"
    T::FT = T_25()
    "Old Temperature `[K]`, if not T, run leaf_temperature_dependence!"
    T_old::FT = 0

    # Photosynthesis system
    "Rate constant for thermal dissipation"
    Kd::FT = 0.85
    "Rate constant for fluorescence (const)"
    Kf::FT = 0.05
    "Reversible NPQ rate constant (initially zero)"
    Kr::FT = 0
    "Sustained NPQ rate constant (for seasonal changes, default is zero)"
    Ks::FT = 0
    "Rate constant for photochemistry (all reaction centers open)"
    Kp::FT = 4
    "Maximal Kp"
    Kp_max::FT = 4
    "max PSII yield (Kr=0, all RC open)"
    maxPSII::FT = Kp / (Kp + Kf + Kd)
    "Fraction of absorbed light used by PSII ETR"
    PSII_frac::FT = 0.5

    # CO₂ pressures
    "Leaf internal CO₂ partial pressure `[Pa]`"
    p_i::FT = 10
    "Leaf surface CO₂ partial pressure `[Pa]`"
    p_s::FT = 40
    "Saturation H₂O vapor pressure `[Pa]`"
    p_sat::FT = saturation_vapor_pressure(T)
    "Leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_bc::FT = 3 / 1.35
    "Leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_lc::FT = 0.01

    # Photosynthesis related
    "RubisCO limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Ac::FT = 0
    "Light limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Aj::FT = 0
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Ag::FT = 0
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    An::FT = 0
    "Product limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    Ap::FT = 0
    "Electron transport `[μmol m⁻² s⁻¹]`"
    J::FT = 0
    "Potential Electron Transport Rate `[μmol m⁻² s⁻¹]`"
    J_pot::FT = 0
    "Maximal electron transport rate `[μmol m⁻² s⁻¹]`"
    Jmax::FT = 120
    "Maximal electron transport rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    Jmax25::FT = 120
    "RubisCO coefficient Kc `[Pa]`"
    Kc::FT = 0
    "RubisCO coefficient Ko `[Pa]`"
    Ko::FT = 0
    "PEP coefficient Ko `[Pa]`"
    Kpep::FT = 0
    "Michaelis-Menten's coefficient `[Pa]`"
    Km::FT = 0
    "Respiration rate `[μmol m⁻² s⁻¹]`"
    Rd::FT = 1
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    Rd25::FT = 1
    "Maximal carboxylation rate `[μmol m⁻² s⁻¹]`"
    Vcmax::FT = 60
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    Vcmax25::FT = 60
    "Maximal PEP carboxylation rate `[μmol m⁻² s⁻¹]`"
    Vpmax::FT = 120
    "Maximal PEP carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    Vpmax25::FT = 120
    "CO₂ compensation point with the absence of Rd `[Pa]`"
    Γ_star::FT = 0

    # Cytochrome related
    "Total concentration of Cytochrome b₆f `[μmol m⁻²]`"
    C_b₆f::FT = 350 / 300
    "Maximal turnover rate of Cytochrome b₆f `[e⁻ s⁻¹]`"
    k_q::FT = 300
    "Maximal Cytochrome b₆f activity `[μmol e⁻ m⁻² s⁻¹]`"
    Vqmax::FT = C_b₆f * k_q
    "Rate constant of consititutive heat loss from the antennae `[s⁻¹]`"
    K_D1::FT = 5.5e8
    "rate constant of fluorescence `[s⁻¹]`"
    K_F1::FT = 5e7
    "Rate constant of regulated heat loss for PS I `[s⁻¹]`"
    K_N1::FT = 1.45e10
    "Rate constant of photochemistry for PS I `[s⁻¹]`"
    K_P1::FT = 1.45e10
    "Rate constant of photochemistry for PS II `[s⁻¹]`"
    K_P2::FT = 4.5e9
    "Rate constant of excitation sharing for PS II `[s⁻¹]`"
    K_U2::FT = 0
    "PPFD absorbed by PS I per incident PPFD"
    α_1::FT = 0.41
    "PPFD absorbed by PS II per incident PPFD"
    α_2::FT = 0.44
    "Weighting factor for PS I"
    ϵ_1::FT = 0
    "Weighting factor for PS II"
    ϵ_2::FT = 1
    "Maximal PS I photochemical yield"
    φ_P1_max::FT = K_P1 / (K_P1 + K_D1 + K_F1)
    "Coupling efficiency of cyclic electron flow `[mol ATP mol⁻¹ e⁻]`"
    n_C::FT = 1
    "Coupling efficiency of linear electron flow `[mol ATP mol⁻¹ e⁻]`"
    n_L::FT = 0.75
    "ratio between J_P700 and J_P680"
    η::FT = 0

    # related to fluorescence
    "PS II electron transport rate `[μmol e⁻ m⁻² s⁻¹]`"
    J_P680_a::FT = 0
    "Rubisco limited PS II electron transport rate `[μmol e⁻ m⁻² s⁻¹]`"
    J_P680_c::FT = 0
    "Light limited PS II electron transport rate `[μmol e⁻ m⁻² s⁻¹]`"
    J_P680_j::FT = 0
    "Product limited PS II electron transport rate `[μmol e⁻ m⁻² s⁻¹]`"
    J_P680_p::FT = 0
    "PS I electron transport rate `[μmol e⁻ m⁻² s⁻¹]`"
    J_P700_a::FT = 0
    "Rubisco limited PS I electron transport rate `[μmol e⁻ m⁻² s⁻¹]`"
    J_P700_c::FT = 0
    "Light limited PS I electron transport rate `[μmol e⁻ m⁻² s⁻¹]`"
    J_P700_j::FT = 0
    "Product limited PS I electron transport rate `[μmol e⁻ m⁻² s⁻¹]`"
    J_P700_p::FT = 0

    # Well watered condition values (to use with β function over PS)
    "Well watered maximal electron transport rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    Jmax25WW::FT = Jmax25
    "Well watered respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    Rd25WW::FT = Rd25
    "Well watered maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    Vcmax25WW::FT = Vcmax25
    "Well watered maximal PEP carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    Vpmax25WW::FT = Vpmax25

    # Fluorescence related
    "Total efficiency, incl. photorespiration `[mol CO₂ mol⁻¹ e-]`"
    e2c::FT = 1 / 6
    "dark adapted yield (`Kp=0`)"
    Fm::FT = 0
    "light adapted yield (`Kp=0`)"
    Fm′::FT = 0
    "dark-adapted fluorescence yield (`Kp=max`)"
    Fo::FT = 0
    "light-adapted fluorescence yield in the dark (`Kp=max`)"
    Fo′::FT = 0
    "Actual electron transport rate `[μmol m⁻² s⁻¹]`"
    Ja::FT = 0
    "Non-Photochemical quenching "
    NPQ::FT = 0
    "Photochemical quenching"
    qQ::FT = 0
    "energy quenching"
    qE::FT = 0
    "PSII yield"
    φ::FT = 0
    "Steady-state (light-adapted) yield (aka Fs)"
    φs::FT = 0

    # Environment related
    "Absorbed photosynthetic active radiation `[μmol m⁻² s⁻¹]`"
    APAR::FT = 100
end








###############################################################################
#
# Photosynthesis modes
#
###############################################################################
"""
    abstract type AbstractCalculationMode

Hierarchy of AbstractCalculationMode
- [`GCO₂Mode`](@ref)
- [`PCO₂Mode`](@ref)
"""
abstract type AbstractCalculationMode end




"""
    struct GCO₂Mode <: AbstractCalculationMode

Calculate leaf photosynthesis using leaf internal CO₂ partial pressure

---
Examples
```julia
# define CO₂ conductance mode
mode = GCO₂Mode();
```
"""
struct GCO₂Mode <: AbstractCalculationMode end




"""
    struct PCO₂Mode <: AbstractCalculationMode

Calculate leaf photosynthesis using total leaf diffusive conductance to CO₂

---
Examples
```julia
# define CO₂ pressure mode
mode = PCO₂Mode();
```
"""
struct PCO₂Mode <: AbstractCalculationMode end








###############################################################################
#
# Leaf fluorescence-related parameter set
#
###############################################################################
"""
    abstract type AbstractFluoModelParaSet{FT}

Hierarchy of the `AbstractFluoModelParaSet`:
- [`FluoParaSet`](@ref)
"""
abstract type AbstractFluoModelParaSet{FT} end




"""
    mutable struct CytoFluoParaSet{FT}

A `AbstractFluoModelParaSet` type paramter set using Johnson-Berry
    photosynthesis model.

---
Examples
```julia
# define Cytochrome model mode (no input)
fluo = CytoFluoParaSet();
```
"""
struct CytoFluoParaSet{FT<:AbstractFloat} <: AbstractFluoModelParaSet{FT} end




"""
    mutable struct FluoParaSet{FT}

A `AbstractFluoModelParaSet` type paramter set.

# Fields
$(TYPEDFIELDS)

---
Examples
```julia
# define classic model mode (3 inputs)
fluo = FluoParaSet(5.01, 1.93, 10.0);
```
"""
struct FluoParaSet{FT<:AbstractFloat} <: AbstractFluoModelParaSet{FT}
    "Fluorescence model coefficient"
    Kr1::FT
    "Fluorescence model coefficient"
    Kr2::FT
    "Fluorescence model coefficient"
    Kr3::FT
end







###############################################################################
#
# Photosynthesis model parameter set -- temperature dependencies and etc
#
###############################################################################
"""
    abstract type AbstractPhotoModelParaSet{FT}

Hierarchy of the `AbstractPhotoModelParaSet`:
- [`C3Cytochrome`](@ref)
- [`C3ParaSet`](@ref)
- [`C4ParaSet`](@ref)
"""
abstract type AbstractPhotoModelParaSet{FT} end




"""
    mutable struct C3Cytochrome{FT}

Parameter sets for C3 photosynthesis with Cytochrome activity.

# Fields
$(TYPEDFIELDS)
"""
mutable struct C3Cytochrome{FT<:AbstractFloat} <: AbstractPhotoModelParaSet{FT}
    "Coefficient 4.0/4.5 for NADPH/ATP requirement stochiometry, respectively"
    Eff_1::FT
    "Coefficient 8.0/10.5 for NADPH/ATP requirement stochiometry, respectively"
    Eff_2::FT
end




"""
    mutable struct C3Paraset{FT}

Parameter sets for C3 photosynthesis.

# Fields
$(TYPEDFIELDS)
"""
mutable struct C3ParaSet{FT<:AbstractFloat} <: AbstractPhotoModelParaSet{FT}
    "Jmax temperature dependency"
    JT::AbstractTemperatureDependency{FT}
    "Kc temperature dependency"
    KcT::AbstractTemperatureDependency{FT}
    "Ko temperature dependency"
    KoT::AbstractTemperatureDependency{FT}
    "Respiration temperature dependency"
    ReT::AbstractTemperatureDependency{FT}
    "Vcmax temperature dependency"
    VcT::AbstractTemperatureDependency{FT}
    "Γ_star temperature dependency"
    ΓsT::AbstractTemperatureDependency{FT}
    "Fluorescence model"
    Flu::AbstractFluoModelParaSet{FT}
    "Coefficient 4.0/4.5 for NADPH/ATP requirement stochiometry, respectively"
    Eff_1::FT
    "Coefficient 8.0/10.5 for NADPH/ATP requirement stochiometry, respectively"
    Eff_2::FT
end




"""
    mutable struct C4ParaSet{FT}

Parameter sets for C3 photosynthesis.

# Fields
$(TYPEDFIELDS)
"""
mutable struct C4ParaSet{FT<:AbstractFloat} <: AbstractPhotoModelParaSet{FT}
    "Kpep temperature dependency"
    KpT::AbstractTemperatureDependency{FT}
    "Respiration temperature dependency"
    ReT::AbstractTemperatureDependency{FT}
    "Vcmax temperature dependency"
    VcT::AbstractTemperatureDependency{FT}
    "Vpmax temperature dependency"
    VpT::AbstractTemperatureDependency{FT}
    "Fluorescence model"
    Flu::AbstractFluoModelParaSet{FT}
end
