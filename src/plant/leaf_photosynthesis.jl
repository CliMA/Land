#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2021-Nov-29: add abstract photosynthesis system type
#     2022-Jan-14: rename to photosynthesis model
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of `AbstractPhotosynthesisModel`:
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
#     2022-Feb-07: add more fields to use with Photosynthesis v0.3.1
#     2022-Feb-07: remove j_p680 and j_p700 series variables
#     2022-Feb-11: split COLIMIT to COLIMIT_CJ and COLIMIT_IP (minor breaking)
#     2022-Mar-01: add two more fields: TD_ΗC and TD_ΗL
#     2022-Mar-01: move η_c and η_l from reaction center to photosynthesis model
#     2022-Mar-09: add v_cmax25_ww to use with StomataModels.jl
#     2022-Jun-13: use Union instead of Abstract... for type definition
#     2022-Jul-18: remove v_cmax25_ww (use β instead)
#     2022-Jul-20: rename TD_ΗC and TD_ΗL to lower greek TD_ηC and TD_ηL
#     2022-Jul-20: move a_c, a_j, and etc as cache variables, and rename them to _a_c, _a_j, and etc
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores C3 Cytochrome photosynthesis system information

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C3CytochromeModel{FT<:AbstractFloat} <: AbstractPhotosynthesisModel{FT}
    # General model information
    "Coefficient 4.0/4.5 for NADPH/ATP requirement stochiometry, respectively"
    EFF_1::FT = 4
    "Coefficient 8.0/10.5 for NADPH/ATP requirement stochiometry, respectively"
    EFF_2::FT = 8

    # Colimitation methods
    "[`AbstractColimit`](@ref) type colimitation method for Ac and Aj => Ai"
    COLIMIT_CJ::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}} = MinimumColimit{FT}()
    "[`AbstractColimit`](@ref) type colimitation method for Ai and Ap => Ag"
    COLIMIT_IP::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}} = MinimumColimit{FT}()
    "[`AbstractColimit`](@ref) type colimitation method for J"
    COLIMIT_J::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}} = SerialColimit{FT}()

    # Temperature dependency structures
    "[`AbstractTemperatureDependency`](@ref) type Kc temperature dependency"
    TD_KC::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = KcTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Ko temperature dependency"
    TD_KO::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = KoTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Kq temperature dependency"
    TD_KQ::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = KqTDJohnson(FT)
    "[`AbstractTemperatureDependency`](@ref) type respiration temperature dependency"
    TD_R::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = RespirationTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Vcmax temperature dependency"
    TD_VCMAX::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = VcmaxTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Γ* temperature dependency"
    TD_Γ::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = ΓStarTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type η_C temperature dependency"
    TD_ηC::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = ηCTDJohnson(FT)
    "[`AbstractTemperatureDependency`](@ref) type η_L temperature dependency"
    TD_ηL::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = ηLTDJohnson(FT)

    # Prognostic variables
    "Total concentration of Cytochrome b₆f `[μmol m⁻²]`"
    b₆f::FT = 350 / 300
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_d25::FT = 0.75
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25::FT = 50

    # Diagnostic variables
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_gross::FT = 0
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_net::FT = 0

    # Cache variables
    "RubisCO limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    _a_c::FT = 0
    "Light limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    _a_j::FT = 0
    "Product limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    _a_p::FT = 0
    "Electron to CO₂ coefficient"
    _e_to_c::FT = 0
    "Potential Electron Transport Rate `[μmol e⁻ m⁻² s⁻¹]`"
    _j_pot::FT = 0
    "PSI electron transport rate after colimitation"
    _j_psi::FT = 0
    "RubisCO coefficient Kc `[Pa]`"
    _k_c::FT = 0
    "Michaelis-Menten's coefficient `[Pa]`"
    _k_m::FT = 0
    "RubisCO coefficient Ko `[Pa]`"
    _k_o::FT = 0
    "Maximal turnover rate of Cytochrome b₆f `[e⁻ s⁻¹]`"
    _k_q::FT = 0
    "Respiration rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    _r_d::FT = 0
    "Last leaf temperature. If different from leaf t, then make temperature correction"
    _t::FT = 0
    "Maximal carboxylation rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    _v_cmax::FT = 0
    "Maximal Cytochrome b₆f activity `[μmol e⁻ m⁻² s⁻¹]`"
    _v_qmax::FT = 0
    "ratio between J_P700 and J_P680"
    _η::FT = 0
    "Coupling efficiency of cyclic electron flow `[mol ATP mol⁻¹ e⁻]`"
    _η_c::FT = 0
    "Coupling efficiency of linear electron flow `[mol ATP mol⁻¹ e⁻]`"
    _η_l::FT = 0
    "CO₂ compensation point with the absence of Rd `[Pa]`"
    _γ_star::FT = 0
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2021-Nov-11: add C3VJPModel structure for classic C₃ photosynthesis system
#     2022-Jan-14: add temperature dependency into the structure
#     2022-Jan-14: rename to photosynthesis model
#     2022-Jan-14: add colimitation and e_to_c
#     2022-Feb-11: split COLIMIT to COLIMIT_CJ, COLIMIT_IP, and COLIMIT_J (minor breaking)
#     2022-Mar-09: add v_cmax25_ww to use with StomataModels.jl
#     2022-Jun-13: use Union instead of Abstract... for type definition
#     2022-Jul-18: remove v_cmax25_ww (use β instead)
#     2022-Jul-20: move a_c, a_j, and etc as cache variables, and rename them to _a_c, _a_j, and etc
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores C3 photosynthesis system information

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C3VJPModel{FT<:AbstractFloat} <: AbstractPhotosynthesisModel{FT}
    # General model information
    "Coefficient 4.0/4.5 for NADPH/ATP requirement stochiometry, respectively"
    EFF_1::FT = 4
    "Coefficient 8.0/10.5 for NADPH/ATP requirement stochiometry, respectively"
    EFF_2::FT = 8

    # Colimitation methods
    "[`AbstractColimit`](@ref) type colimitation method for Ac and Aj => Ai"
    COLIMIT_CJ::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}} = MinimumColimit{FT}()
    "[`AbstractColimit`](@ref) type colimitation method for Ai and Ap => Ag"
    COLIMIT_IP::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}} = MinimumColimit{FT}()
    "[`AbstractColimit`](@ref) type colimitation method for J"
    COLIMIT_J::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}} = MinimumColimit{FT}()

    # Temperature dependency structures
    "[`AbstractTemperatureDependency`](@ref) type Jmax temperature dependency"
    TD_JMAX::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = JmaxTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Kc temperature dependency"
    TD_KC::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = KcTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Ko temperature dependency"
    TD_KO::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = KoTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type respiration temperature dependency"
    TD_R::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = RespirationTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Vcmax temperature dependency"
    TD_VCMAX::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = VcmaxTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Γ* temperature dependency"
    TD_Γ::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = ΓStarTDCLM(FT)

    # Prognostic variables
    "Maximal electron transport rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    j_max25::FT = 83.5
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_d25::FT = 0.75
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25::FT = 50

    # Diagnostic variables
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_gross::FT = 0
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_net::FT = 0

    # Cache variables
    "RubisCO limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    _a_c::FT = 0
    "Light limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    _a_j::FT = 0
    "Product limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    _a_p::FT = 0
    "Electron to CO₂ coefficient"
    _e_to_c::FT = 0
    "Electron transport `[μmol m⁻² s⁻¹]`"
    _j::FT = 0
    "Maximal electron transport rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    _j_max::FT = 0
    "Potential Electron Transport Rate `[μmol m⁻² s⁻¹]`"
    _j_pot::FT = 0
    "RubisCO coefficient Kc `[Pa]`"
    _k_c::FT = 0
    "Michaelis-Menten's coefficient `[Pa]`"
    _k_m::FT = 0
    "RubisCO coefficient Ko `[Pa]`"
    _k_o::FT = 0
    "Respiration rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    _r_d::FT = 0
    "Last leaf temperature. If different from leaf t, then make temperature correction"
    _t::FT = 0
    "Maximal carboxylation rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    _v_cmax::FT = 0
    "CO₂ compensation point with the absence of Rd `[Pa]`"
    _γ_star::FT = 0
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-14: add C4VJPModel structure for classic C₄ photosynthesis system
#     2022-Feb-11: remove j from the struct
#     2022-Feb-11: split COLIMIT to COLIMIT_CJ and COLIMIT_IP (minor breaking)
#     2022-Mar-09: add v_cmax25_ww to use with StomataModels.jl
#     2022-Jun-13: use Union instead of Abstract... for type definition
#     2022-Jul-18: remove v_cmax25_ww (use β instead)
#     2022-Jul-20: move a_c, a_j, and etc as cache variables, and rename them to _a_c, _a_j, and etc
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
Base.@kwdef mutable struct C4VJPModel{FT<:AbstractFloat} <: AbstractPhotosynthesisModel{FT}
    # Colimitation methods
    "[`AbstractColimit`](@ref) type colimitation method for Ac and Aj => Ai"
    COLIMIT_CJ::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}} = MinimumColimit{FT}()
    "[`AbstractColimit`](@ref) type colimitation method for Ai and Ap => Ag"
    COLIMIT_IP::Union{MinimumColimit{FT}, QuadraticColimit{FT}, SerialColimit{FT}} = MinimumColimit{FT}()

    # Temperature dependency structures
    "[`AbstractTemperatureDependency`](@ref) type Kpep temperature dependency"
    TD_KPEP::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = KpepTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type  respiration temperature dependency"
    TD_R::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = RespirationTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Vcmax temperature dependency"
    TD_VCMAX::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = VcmaxTDCLM(FT)
    "[`AbstractTemperatureDependency`](@ref) type Vpmax temperature dependency"
    TD_VPMAX::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}} = VpmaxTDBoyd(FT)

    # Prognostic variables
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_d25::FT = 0.75
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25::FT = 50
    "Maximal PEP carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_pmax25::FT = 50

    # Diagnostic variables
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_gross::FT = 0
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_net::FT = 0

    # Cache variables
    "RubisCO limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    _a_c::FT = 0
    "Light limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    _a_j::FT = 0
    "Product limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    _a_p::FT = 0
    "Electron to CO₂ coefficient"
    _e_to_c::FT = 0
    "Potential Electron Transport Rate `[μmol m⁻² s⁻¹]`"
    _j_pot::FT = 0
    "PEP coefficient Kpep `[Pa]`"
    _k_pep::FT = 0
    "Respiration rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    _r_d::FT = 0
    "Last leaf temperature. If different from leaf t, then make temperature correction"
    _t::FT = 0
    "Maximal carboxylation rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    _v_cmax::FT = 0
    "Maximal PEP carboxylation rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    _v_pmax::FT = 0
end


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jan-14: add abstract mode type
#     2022-Jan-14: add Hierarchy description
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractPhotosynthesisMode:
- [`GCO₂Mode`](@ref)
- [`PCO₂Mode`](@ref)

"""
abstract type AbstractPhotosynthesisMode end


""" An empty structure to signal the function to calculate photosynthetic rates based on leaf diffusive conductance to CO₂ """
struct GCO₂Mode <: AbstractPhotosynthesisMode end


""" An empty structure to signal the function to calculate photosynthetic rates based on CO₂ partial pressure """
struct PCO₂Mode <: AbstractPhotosynthesisMode end
