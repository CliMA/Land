"""
$(TYPEDEF)

Hierachy of AbstractSoilVC:
- [`C3VJPModel`](@ref)
- [`C4VJPModel`](@ref)
"""
abstract type AbstractPhotosynthesisModel{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2021-Nov-11: add C3VJPModel structure for classic C₃ photosynthesis system
#     2022-Jan-14: add temperature dependency into the structure
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
    "Colimitation method"
    COLIMIT::AbstractColimit{FT}
    "Coefficient 4.0/4.5 for NADPH/ATP requirement stochiometry, respectively"
    EFF_1::FT
    "Coefficient 8.0/10.5 for NADPH/ATP requirement stochiometry, respectively"
    EFF_2::FT
    "Jmax temperature dependency"
    TD_JMAX::AbstractTemperatureDependency{FT}
    "Kc temperature dependency"
    TD_KC::AbstractTemperatureDependency{FT}
    "Ko temperature dependency"
    TD_KO::AbstractTemperatureDependency{FT}
    "Respiration temperature dependency"
    TD_R::AbstractTemperatureDependency{FT}
    "Vcmax temperature dependency"
    TD_VCMAX::AbstractTemperatureDependency{FT}
    "Γ* temperature dependency"
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
#     2021-Jan-14: add C4VJPModel structure for classic C₄ photosynthesis system
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
    "Colimitation method"
    COLIMIT::AbstractColimit{FT}
    "Kpep temperature dependency"
    TD_KPEP::AbstractTemperatureDependency{FT}
    "Respiration temperature dependency"
    TD_R::AbstractTemperatureDependency{FT}
    "Vcmax temperature dependency"
    TD_VCMAX::AbstractTemperatureDependency{FT}
    "Vpmax temperature dependency"
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


"""
C4VJPModel{FT}(; v_cmax25::Number = 50, j_max25::Number = 83.5, r_d25::Number = 0.75) where {FT<:AbstractFloat}

Constructor for [`C4VJPModel`](@ref), given
- `v_cmax25` Maximal carboxylation rate at 298.15 K
- `v_pmax25` Maximal PEP carboxylation rate at 298.15 K
- `r_d25` Respiration rate at 298.15 K

---
# Examples
```julia
c4 = C4VJPModel{Float64}();
c4 = C4VJPModel{Float64}(v_cmax25 = 30, v_pmax = 40, r_d25 = 1);
```
"""
C4VJPModel{FT}(; v_cmax25::Number = 50, v_pmax25::Number = 50, r_d25::Number = 0.75) where {FT<:AbstractFloat} = (
    return C4VJPModel{FT}(MinimumColimit{FT}(), KpepTDCLM(FT), RespirationTDCLM(FT), VcmaxTDCLM(FT), VpmaxTDBoyd(FT),
                          r_d25, v_cmax25, v_pmax25, 0, 0, 0, -r_d25, 0, 0, 0, 0, 0, r_d25, v_cmax25, v_pmax25)
);


"""
$(TYPEDEF)

Hierachy of AbstractSoilVC:
- [`GCO₂Mode`](@ref)
- [`PCO₂Mode`](@ref)
"""
abstract type AbstractPhotosynthesisMode end


"""
$(TYPEDEF)

An empty structure to signal the function to calculate photosynthetic rates based on leaf diffusive conductance to CO₂.
"""
struct GCO₂Mode <: AbstractPhotosynthesisMode end


"""
$(TYPEDEF)

An empty structure to signal the function to calculate photosynthetic rates based on CO₂ partial pressure.
"""
struct PCO₂Mode <: AbstractPhotosynthesisMode end
