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
#     2022-Feb-07: add more fields to use with Photosynthesis v0.3.1
#     2022-Feb-07: remove j_p680 and j_p700 series variables
#     2022-Feb-11: split COLIMIT to COLIMIT_CJ and COLIMIT_IP (minor breaking)
#     2022-Mar-01: add two more fields: TD_ΗC and TD_ΗL
#     2022-Mar-01: move η_c and η_l from reaction center to photosynthesis model
#     2022-Mar-01: move k_q from prognostic to dignostic section
#     2022-Mar-09: add v_cmax25_ww to use with StomataModels.jl
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
    "[`AbstractColimit`](@ref) type colimitation method for Ac and Aj => Ai"
    COLIMIT_CJ::AbstractColimit{FT}
    "[`AbstractColimit`](@ref) type colimitation method for Ai and Ap => Ag"
    COLIMIT_IP::AbstractColimit{FT}
    "[`AbstractColimit`](@ref) type colimitation method for J"
    COLIMIT_J::AbstractColimit{FT}
    "Coefficient 4.0/4.5 for NADPH/ATP requirement stochiometry, respectively"
    EFF_1::FT
    "Coefficient 8.0/10.5 for NADPH/ATP requirement stochiometry, respectively"
    EFF_2::FT
    "[`AbstractTemperatureDependency`](@ref) type Kc temperature dependency"
    TD_KC::AbstractTemperatureDependency{FT}
    "[`AbstractTemperatureDependency`](@ref) type Ko temperature dependency"
    TD_KO::AbstractTemperatureDependency{FT}
    "[`AbstractTemperatureDependency`](@ref) type Kq temperature dependency"
    TD_KQ::AbstractTemperatureDependency{FT}
    "[`AbstractTemperatureDependency`](@ref) type respiration temperature dependency"
    TD_R::AbstractTemperatureDependency{FT}
    "[`AbstractTemperatureDependency`](@ref) type Vcmax temperature dependency"
    TD_VCMAX::AbstractTemperatureDependency{FT}
    "[`AbstractTemperatureDependency`](@ref) type Γ* temperature dependency"
    TD_Γ::AbstractTemperatureDependency{FT}
    "[`AbstractTemperatureDependency`](@ref) type Η_C temperature dependency"
    TD_ΗC::AbstractTemperatureDependency{FT}
    "[`AbstractTemperatureDependency`](@ref) type Η_L temperature dependency"
    TD_ΗL::AbstractTemperatureDependency{FT}

    # prognostic variables that change with time
    "Total concentration of Cytochrome b₆f `[μmol m⁻²]`"
    b₆f::FT
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_d25::FT
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25::FT
    "Well watered maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25_ww::FT

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
    "Potential Electron Transport Rate `[μmol e⁻ m⁻² s⁻¹]`"
    j_pot::FT
    "PSI electron transport rate after colimitation"
    j_psi::FT
    "RubisCO coefficient Kc `[Pa]`"
    k_c::FT
    "Michaelis-Menten's coefficient `[Pa]`"
    k_m::FT
    "RubisCO coefficient Ko `[Pa]`"
    k_o::FT
    "Maximal turnover rate of Cytochrome b₆f `[e⁻ s⁻¹]`"
    k_q::FT
    "Respiration rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    r_d::FT
    "Maximal carboxylation rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    v_cmax::FT
    "Maximal Cytochrome b₆f activity `[μmol e⁻ m⁻² s⁻¹]`"
    v_qmax::FT
    "ratio between J_P700 and J_P680"
    η::FT
    "Coupling efficiency of cyclic electron flow `[mol ATP mol⁻¹ e⁻]`"
    η_c::FT
    "Coupling efficiency of linear electron flow `[mol ATP mol⁻¹ e⁻]`"
    η_l::FT
    "CO₂ compensation point with the absence of Rd `[Pa]`"
    γ_star::FT
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2021-Nov-18: add constructor
#     2022-Jan-25: add documentation
#     2022-Feb-07: add more fields into the constructors
#     2022-Feb-11: split COLIMIT to COLIMIT_CJ and COLIMIT_IP (minor breaking)
#     2022-Feb-11: add colimit option in constructor to enable quick deployment of quadratic colimitation
#     2022-Mar-01: add two more fields: TD_ΗC and TD_ΗL
#     2022-Mar-01: move η_c and η_l from reaction center to photosynthesis model
#     2022-Mar-01: move k_q from prognostic to dignostic section
#     2022-Mar-09: add v_cmax25_ww to use with StomataModels.jl
#
#######################################################################################################################################################################################################
"""

    C3CytochromeModel{FT}(; v_cmax25::Number = 50, r_d25::Number = 0.75, colimit::Bool = false) where {FT<:AbstractFloat}

Constructor for `C3CytochromeModel`, given
- `v_cmax25` Maximal carboxylation rate at 298.15 K
- `r_d25` Respiration rate at 298.15 K
- `colimit` If true, use quadratic colimitations for a_c, a_j, and a_p

---
# Examples
```julia
cy = C3CytochromeModel{Float64}();
cy = C3CytochromeModel{Float64}(v_cmax25 = 30, r_d25 = 1, colimit = true);
```
"""
C3CytochromeModel{FT}(; v_cmax25::Number = 50, r_d25::Number = 0.75, colimit::Bool = false) where {FT<:AbstractFloat} = (
    if colimit
        _colim_cj = QuadraticColimit{FT}(0.98);
        _colim_ip = QuadraticColimit{FT}(0.95);
    else
        _colim_cj = MinimumColimit{FT}();
        _colim_ip = MinimumColimit{FT}();
    end;

    return C3CytochromeModel{FT}(
                _colim_cj,              # COLIMIT_CJ
                _colim_ip,              # COLIMIT_IP
                SerialColimit{FT}(),    # COLIMIT_J
                4,                      # EFF_1
                8,                      # EFF_2
                KcTDCLM(FT),            # TD_KC
                KoTDCLM(FT),            # TD_KO
                KqTDJohnson(FT),        # TD_KQ
                RespirationTDCLM(FT),   # TD_R
                VcmaxTDCLM(FT),         # TD_VCMAX
                ΓStarTDCLM(FT),         # TD_Γ
                ΗCTDJohnson(FT),        # TD_ΗC
                ΗLTDJohnson(FT),        # TD_ΗL
                350 / 300,              # b₆f
                r_d25,                  # r_d25
                v_cmax25,               # v_cmax25,
                v_cmax25,               # v_cmax25_ww
                0,                      # a_c
                0,                      # a_gross
                0,                      # a_j
                -r_d25,                 # a_net
                0,                      # a_p
                0,                      # e_to_c
                0,                      # j_pot
                0,                      # j_psi
                0,                      # k_c
                0,                      # k_m
                0,                      # k_o
                0,                      # k_q
                r_d25,                  # r_d
                v_cmax25,               # v_cmax
                0,                      # v_qmax
                0,                      # η
                0,                      # η_c
                0,                      # η_l
                0)                      # γ_star
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
#     2022-Feb-11: split COLIMIT to COLIMIT_CJ, COLIMIT_IP, and COLIMIT_J (minor breaking)
#     2022-Mar-09: add v_cmax25_ww to use with StomataModels.jl
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
    "[`AbstractColimit`](@ref) type colimitation method for Ac and Aj => Ai"
    COLIMIT_CJ::AbstractColimit{FT}
    "[`AbstractColimit`](@ref) type colimitation method for Ai and Ap => Ag"
    COLIMIT_IP::AbstractColimit{FT}
    "[`AbstractColimit`](@ref) type colimitation method for J"
    COLIMIT_J::AbstractColimit{FT}
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
    "Well watered maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25_ww::FT

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
#     2022-Feb-11: split COLIMIT to COLIMIT_CJ, COLIMIT_IP, and COLIMIT_J (minor breaking)
#     2022-Feb-11: add colimit option in constructor to enable quick deployment of quadratic colimitation
#     2022-Mar-09: add v_cmax25_ww to use with StomataModels.jl
#
#######################################################################################################################################################################################################
"""

    C3VJPModel{FT}(; v_cmax25::Number = 50, j_max25::Number = 83.5, r_d25::Number = 0.75, colimit::Bool = false) where {FT<:AbstractFloat}

Constructor for [`C3VJPModel`](@ref), given
- `v_cmax25` Maximal carboxylation rate at 298.15 K
- `j_max25` Maximal electron transport rate at 298.15 K
- `r_d25` Respiration rate at 298.15 K
- `colimit` If true, use quadratic colimitations for j and a_c, a_j, and a_p

---
# Examples
```julia
c3 = C3VJPModel{Float64}();
c3 = C3VJPModel{Float64}(v_cmax25 = 30, j_max25 = 50, r_d25 = 1, colimit = true);
```
"""
C3VJPModel{FT}(; v_cmax25::Number = 50, j_max25::Number = 83.5, r_d25::Number = 0.75, colimit::Bool = false) where {FT<:AbstractFloat} = (
    if colimit
        _colim_cj = QuadraticColimit{FT}(0.98);
        _colim_ip = QuadraticColimit{FT}(0.95);
        _colim_j  = QuadraticColimit{FT}(0.7);
    else
        _colim_cj = MinimumColimit{FT}();
        _colim_ip = MinimumColimit{FT}();
        _colim_j  = MinimumColimit{FT}();
    end;

    return C3VJPModel{FT}(
                _colim_cj,              # COLIMIT_CJ
                _colim_ip,              # COLIMIT_IP
                _colim_j,               # COLIMIT_J
                4,                      # EFF_1
                8,                      # EFF_2
                JmaxTDCLM(FT),          # TD_JMAX
                KcTDCLM(FT),            # TD_KC
                KoTDCLM(FT),            # TD_KO
                RespirationTDCLM(FT),   # TD_R
                VcmaxTDCLM(FT),         # TD_VCMAX
                ΓStarTDCLM(FT),         # TD_Γ
                j_max25,                # j_max25
                r_d25,                  # r_d25
                v_cmax25,               # v_cmax25
                v_cmax25,               # v_cmax25_ww
                0,                      # a_c
                0,                      # a_gross
                0,                      # a_j
                -r_d25,                 # a_net
                0,                      # a_p
                0,                      # e_to_c
                0,                      # j
                j_max25,                # j_max
                0,                      # j_pot
                0,                      # k_c
                0,                      # k_m
                0,                      # k_o
                r_d25,                  # r_d
                v_cmax25,               # v_cmax
                0)                      # γ_star
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-14: add C4VJPModel structure for classic C₄ photosynthesis system
#     2022-Jan-25: fix documentation
#     2022-Feb-11: remove j from the struct
#     2022-Feb-11: split COLIMIT to COLIMIT_CJ and COLIMIT_IP (minor breaking)
#     2022-Mar-09: add v_cmax25_ww to use with StomataModels.jl
# To do
#     TODO: add Jmax to C4VJPModel and thus JMAX TD in Photosynthesis.jl (not necessary)
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
    "[`AbstractColimit`](@ref) type colimitation method for Ac and Aj => Ai"
    COLIMIT_CJ::AbstractColimit{FT}
    "[`AbstractColimit`](@ref) type colimitation method for Ai and Ap => Ag"
    COLIMIT_IP::AbstractColimit{FT}
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
    "Well watered maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25_ww::FT
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
#     2022-Feb-11: remove j from the struct
#     2022-Feb-11: default e_to_c set to 1/6
#     2022-Feb-11: split COLIMIT to COLIMIT_CJ and COLIMIT_IP (minor breaking)
#     2022-Feb-11: add colimit option in constructor to enable quick deployment of quadratic colimitation
#     2022-Mar-09: add v_cmax25_ww to use with StomataModels.jl
#
#######################################################################################################################################################################################################
"""

    C4VJPModel{FT}(; v_cmax25::Number = 50, v_pmax25::Number = 50, r_d25::Number = 0.75, colimit::Bool = false) where {FT<:AbstractFloat}

Constructor for [`C4VJPModel`](@ref), given
- `v_cmax25` Maximal carboxylation rate at 298.15 K
- `v_pmax25` Maximal PEP carboxylation rate at 298.15 K
- `r_d25` Respiration rate at 298.15 K
- `colimit` If true, use quadratic colimitations for a_c, a_j, and a_p

---
# Examples
```julia
c4 = C4VJPModel{Float64}();
c4 = C4VJPModel{Float64}(v_cmax25 = 30, v_pmax25 = 40, r_d25 = 1, colimit = true);
```
"""
C4VJPModel{FT}(; v_cmax25::Number = 50, v_pmax25::Number = 50, r_d25::Number = 0.75, colimit::Bool = false) where {FT<:AbstractFloat} = (
    if colimit
        _colim_cj = QuadraticColimit{FT}(0.8);
        _colim_ip = QuadraticColimit{FT}(0.95);
    else
        _colim_cj = MinimumColimit{FT}();
        _colim_ip = MinimumColimit{FT}();
    end;

    return C4VJPModel{FT}(
                _colim_cj,              # COLIMIT_CJ
                _colim_ip,              # COLIMIT_IP
                KpepTDCLM(FT),          # TD_KPEP
                RespirationTDCLM(FT),   # TD_R
                VcmaxTDCLM(FT),         # TD_VCMAX
                VpmaxTDBoyd(FT),        # TD_VPMAX
                r_d25,                  # r_d25
                v_cmax25,               # v_cmax25
                v_cmax25,               # v_cmax25_ww
                v_pmax25,               # v_pmax25
                0,                      # a_c
                0,                      # a_gross
                0,                      # a_j
                -r_d25,                 # a_net
                0,                      # a_p
                1/6,                    # e_to_c
                0,                      # j_pot
                0,                      # k_pep
                r_d25,                  # r_d
                v_cmax25,               # v_cmax
                v_pmax25)               # v_pmax
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
