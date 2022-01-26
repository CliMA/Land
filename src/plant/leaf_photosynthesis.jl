#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2021-Nov-29: add abstract photosynthesis system type
#     2022-Jan-14: rename to photosynthesis model
#     2022-Jan-14: add hierachy description
#     2022-Jan-25: fix documentation
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierachy of `AbstractPhotosynthesisModel`:
- [`C3CytochromeModel`](@ref)
- [`C3VJPModel`](@ref)
- [`C4VJPModel`](@ref)
"""
abstract type AbstractPhotosynthesisModel{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-18: add C3CytochromeModel structure for C₃ photosynthesis system
#     2022-Jan-25: fix documentation
# To do
#     TODO: add TD in Photosynthesis.jl
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores C3 Cytochrome photosynthesis system information

# Fields

$(TYPEDFIELDS)

"""
mutable struct C3CytochromeModel{FT<:AbstractFloat} <: AbstractPhotosynthesisModel{FT}
    # parameters that do not change with time
    "[`AbstractColimit`](@ref) type colimitation method"
    COLIMIT::AbstractColimit{FT}
    "Coefficient 4.0/4.5 for NADPH/ATP requirement stochiometry, respectively"
    EFF_1::FT
    "Coefficient 8.0/10.5 for NADPH/ATP requirement stochiometry, respectively"
    EFF_2::FT

    # prognostic variables that change with time
    "Total concentration of Cytochrome b₆f `[μmol m⁻²]`"
    b₆f::FT
    "Maximal turnover rate of Cytochrome b₆f `[e⁻ s⁻¹]`"
    k_q::FT

    # dignostic variables that change with time
    "PS II electron transport rate `[μmol e⁻ m⁻² s⁻¹]`"
    j_p680_a::FT
    "Rubisco limited PS II electron transport rate `[μmol e⁻ m⁻² s⁻¹]`"
    j_p680_c::FT
    "Light limited PS II electron transport rate `[μmol e⁻ m⁻² s⁻¹]`"
    j_p680_j::FT
    "Product limited PS II electron transport rate `[μmol e⁻ m⁻² s⁻¹]`"
    j_p680_p::FT
    "PS I electron transport rate `[μmol e⁻ m⁻² s⁻¹]`"
    j_p700_a::FT
    "Rubisco limited PS I electron transport rate `[μmol e⁻ m⁻² s⁻¹]`"
    j_p700_c::FT
    "Light limited PS I electron transport rate `[μmol e⁻ m⁻² s⁻¹]`"
    j_p700_j::FT
    "Product limited PS I electron transport rate `[μmol e⁻ m⁻² s⁻¹]`"
    j_p700_p::FT
    "Maximal Cytochrome b₆f activity `[μmol e⁻ m⁻² s⁻¹]`"
    v_qmax::FT
    "ratio between J_P700 and J_P680"
    η::FT
    "CO₂ compensation point with the absence of Rd `[Pa]`"
    γ_star::FT
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2021-Nov-18: add constructor
#     2022-Jan-25: add documentation
#
#######################################################################################################################################################################################################
"""

    C3CytochromeModel{FT}(; v_cmax25::Number = 50, r_d25::Number = 0.75) where {FT<:AbstractFloat}

Constructor for `C3CytochromeModel`, given
- `v_cmax25` Maximal carboxylation rate at 298.15 K
- `r_d25` Respiration rate at 298.15 K

---
# Examples
```julia
cy = C3CytochromeModel{Float64}();
cy = C3CytochromeModel{Float64}(v_cmax25 = 30, r_d25 = 1);
```
"""
C3CytochromeModel{FT}(; v_cmax25::Number = 50, r_d25::Number = 0.75) where {FT<:AbstractFloat} = (
    return C3CytochromeModel{FT}(MinimumColimit{FT}(), 4, 8, 350 / 300, 300, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2021-Nov-11: add C3VJPModel structure for classic C₃ photosynthesis system
#     2022-Jan-14: add temperature dependency into the structure
#     2022-Jan-14: rename to photosynthesis model
#     2022-Jan-14: add colimitation and e_to_c
#     2022-Jan-25: fix documentation
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores C3 photosynthesis system information

# Fields

$(TYPEDFIELDS)

"""
mutable struct C3VJPModel{FT<:AbstractFloat} <: AbstractPhotosynthesisModel{FT}
    # parameters that do not change with time
    "[`AbstractColimit`](@ref) type colimitation method"
    COLIMIT::AbstractColimit{FT}
    "Coefficient 4.0/4.5 for NADPH/ATP requirement stochiometry, respectively"
    EFF_1::FT
    "Coefficient 8.0/10.5 for NADPH/ATP requirement stochiometry, respectively"
    EFF_2::FT
    "[`AbstractTemperatureDependency`](@ref) type Jmax temperature dependency"
    TD_JMAX::AbstractTemperatureDependency{FT}
    "[`AbstractTemperatureDependency`](@ref) type Kc temperature dependency"
    TD_KC::AbstractTemperatureDependency{FT}
    "[`AbstractTemperatureDependency`](@ref) type Ko temperature dependency"
    TD_KO::AbstractTemperatureDependency{FT}
    "[`AbstractTemperatureDependency`](@ref) type respiration temperature dependency"
    TD_R::AbstractTemperatureDependency{FT}
    "[`AbstractTemperatureDependency`](@ref) type Vcmax temperature dependency"
    TD_VCMAX::AbstractTemperatureDependency{FT}
    "[`AbstractTemperatureDependency`](@ref) type Γ* temperature dependency"
    TD_Γ::AbstractTemperatureDependency{FT}

    # prognostic variables that change with time
    "Maximal electron transport rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    j_max25::FT
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_d25::FT
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25::FT

    # dignostic variables that change with time
    "RubisCO limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_c::FT
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_gross::FT
    "Light limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_j::FT
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_net::FT
    "Product limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_p::FT
    "Electron to CO₂ coefficient"
    e_to_c::FT
    "Electron transport `[μmol m⁻² s⁻¹]`"
    j::FT
    "Maximal electron transport rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    j_max::FT
    "Potential Electron Transport Rate `[μmol m⁻² s⁻¹]`"
    j_pot::FT
    "RubisCO coefficient Kc `[Pa]`"
    k_c::FT
    "Michaelis-Menten's coefficient `[Pa]`"
    k_m::FT
    "RubisCO coefficient Ko `[Pa]`"
    k_o::FT
    "Respiration rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    r_d::FT
    "Maximal carboxylation rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    v_cmax::FT
    "CO₂ compensation point with the absence of Rd `[Pa]`"
    γ_star::FT
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2021-Nov-11: add constructor
#     2022-Jan-14: update constructor to match structure
#     2022-Jan-25: fix documentation
#
#######################################################################################################################################################################################################
"""

    C3VJPModel{FT}(; v_cmax25::Number = 50, j_max25::Number = 83.5, r_d25::Number = 0.75) where {FT<:AbstractFloat}

Constructor for [`C3VJPModel`](@ref), given
- `v_cmax25` Maximal carboxylation rate at 298.15 K
- `j_max25` Maximal electron transport rate at 298.15 K
- `r_d25` Respiration rate at 298.15 K

---
# Examples
```julia
c3 = C3VJPModel{Float64}();
c3 = C3VJPModel{Float64}(v_cmax25 = 30, j_max25 = 50, r_d25 = 1);
```
"""
C3VJPModel{FT}(; v_cmax25::Number = 50, j_max25::Number = 83.5, r_d25::Number = 0.75) where {FT<:AbstractFloat} = (
    return C3VJPModel{FT}(MinimumColimit{FT}(), 4, 8, JmaxTDCLM(FT), KcTDCLM(FT), KoTDCLM(FT), RespirationTDCLM(FT), VcmaxTDCLM(FT), ΓStarTDCLM(FT),
                          j_max25, r_d25, v_cmax25, 0, 0, 0, -r_d25, 0, 0, 0, j_max25, 0, 0, 0, 0, r_d25, v_cmax25, 0)
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-14: add C4VJPModel structure for classic C₄ photosynthesis system
#     2022-Jan-25: fix documentation
# To do
#     TODO: add Jmax to C4VJPModel and thus JMAX TD in Photosynthesis.jl
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores C4 photosynthesis system information

# Fields

$(TYPEDFIELDS)

"""
mutable struct C4VJPModel{FT<:AbstractFloat} <: AbstractPhotosynthesisModel{FT}
    # parameters that do not change with time
    "[`AbstractColimit`](@ref) type colimitation method"
    COLIMIT::AbstractColimit{FT}
    "[`AbstractTemperatureDependency`](@ref) type Kpep temperature dependency"
    TD_KPEP::AbstractTemperatureDependency{FT}
    "[`AbstractTemperatureDependency`](@ref) type  respiration temperature dependency"
    TD_R::AbstractTemperatureDependency{FT}
    "[`AbstractTemperatureDependency`](@ref) type Vcmax temperature dependency"
    TD_VCMAX::AbstractTemperatureDependency{FT}
    "[`AbstractTemperatureDependency`](@ref) type Vpmax temperature dependency"
    TD_VPMAX::AbstractTemperatureDependency{FT}

    # prognostic variables that change with time
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_d25::FT
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25::FT
    "Maximal PEP carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_pmax25::FT

    # dignostic variables that change with time
    "RubisCO limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_c::FT
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_gross::FT
    "Light limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_j::FT
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_net::FT
    "Product limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_p::FT
    "Electron to CO₂ coefficient"
    e_to_c::FT
    "Electron transport `[μmol m⁻² s⁻¹]`"
    j::FT
    "Potential Electron Transport Rate `[μmol m⁻² s⁻¹]`"
    j_pot::FT
    "PEP coefficient Kpep `[Pa]`"
    k_pep::FT
    "Respiration rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    r_d::FT
    "Maximal carboxylation rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    v_cmax::FT
    "Maximal PEP carboxylation rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    v_pmax::FT
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2021-Nov-18: add constructor
#     2022-Jan-25: fix documentation
#
#######################################################################################################################################################################################################
"""

    C4VJPModel{FT}(; v_cmax25::Number = 50, v_pmax25::Number = 50, r_d25::Number = 0.75) where {FT<:AbstractFloat}

Constructor for [`C4VJPModel`](@ref), given
- `v_cmax25` Maximal carboxylation rate at 298.15 K
- `v_pmax25` Maximal PEP carboxylation rate at 298.15 K
- `r_d25` Respiration rate at 298.15 K

---
# Examples
```julia
c4 = C4VJPModel{Float64}();
c4 = C4VJPModel{Float64}(v_cmax25 = 30, v_pmax25 = 40, r_d25 = 1);
```
"""
C4VJPModel{FT}(; v_cmax25::Number = 50, v_pmax25::Number = 50, r_d25::Number = 0.75) where {FT<:AbstractFloat} = (
    return C4VJPModel{FT}(MinimumColimit{FT}(), KpepTDCLM(FT), RespirationTDCLM(FT), VcmaxTDCLM(FT), VpmaxTDBoyd(FT),
                          r_d25, v_cmax25, v_pmax25, 0, 0, 0, -r_d25, 0, 0, 0, 0, 0, r_d25, v_cmax25, v_pmax25)
);


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jan-14: add abstract mode type
#     2022-Jan-14: add hierachy description
#     2022-Jan-25: fix documentation
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierachy of AbstractSoilVC:
- [`GCO₂Mode`](@ref)
- [`PCO₂Mode`](@ref)
"""
abstract type AbstractPhotosynthesisMode end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-14: add CO₂ conductance mode
#     2022-Jan-25: fix documentation
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

An empty structure to signal the function to calculate photosynthetic rates based on leaf diffusive conductance to CO₂.

---
# Examples
```julia
mode = GCO₂Mode();
```
"""
struct GCO₂Mode <: AbstractPhotosynthesisMode end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-14: add CO₂ partial pressure mode
#     2022-Jan-25: fix documentation
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

An empty structure to signal the function to calculate photosynthetic rates based on CO₂ partial pressure.

---
# Examples
```julia
mode = PCO₂Mode();
```
"""
struct PCO₂Mode <: AbstractPhotosynthesisMode end
