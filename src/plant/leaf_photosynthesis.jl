"""
$(TYPEDEF)

Hierachy of AbstractSoilVC:
- [`C3VJPSystem`](@ref)
"""
abstract type AbstractPhotosynthesisSystem{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2021-Nov-11: add C3VJPSystem structure for classic C₃ photosynthesis system
#
#######################################################################################################################################################################################################
mutable struct C3VJPSystem{FT<:AbstractFloat} <: AbstractPhotosynthesisSystem{FT}
    # parameters that do not change with time
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
    "Coefficient 4.0/4.5 for NADPH/ATP requirement stochiometry, respectively"
    EFF_1::FT
    "Coefficient 8.0/10.5 for NADPH/ATP requirement stochiometry, respectively"
    EFF_2::FT

    # prognostic variables that change with time
    "Maximal electron transport rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    j_max25::FT
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_d25::FT
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_max25::FT

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
    v_max::FT
    "CO₂ compensation point with the absence of Rd `[Pa]`"
    γ_star::FT
end


"""
    C3VJPSystem{FT}(; v_max25::Number = 50, j_max25::Number = 83.5, r_d25::Number = 0.75) where {FT<:AbstractFloat}

Constructor for [`C3VJPSystem`](@ref), given
- `v_max25` Maximal carboxylation rate at 298.15 K
- `j_max25` Maximal electron transport rate at 298.15 K
- `r_d25` Respiration rate at 298.15 K

---
# Examples
```julia
c3 = C3VJPSystem{Float64}();
c3 = C3VJPSystem{Float64}(v_max25 = 30, j_max25 = 50, r_d25 = 1);
```
"""
C3VJPSystem{FT}(; v_max25::Number = 50, j_max25::Number = 83.5, r_d25::Number = 0.75) where {FT<:AbstractFloat} = (
    return C3VJPSystem{FT}(JmaxTDCLM(FT), KcTDCLM(FT), KoTDCLM(FT), RespirationTDCLM(FT), VcmaxTDCLM(FT),ΓStarTDCLM(FT), 4, 8,
                           j_max25, r_d25, v_max25, 0, 0, -r_d25, 0, 0, 0, j_max25, 0, 0, 0, 0, r_d25, v_max25, 0)
);
