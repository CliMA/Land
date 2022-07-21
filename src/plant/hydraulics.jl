#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-May-24: add abstract type for hydraulic organ
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractHydraulicSystem:
- [`LeafHydraulics`](@ref)
- [`RootHydraulics`](@ref)
- [`StemHydraulics`](@ref)

"""
abstract type AbstractHydraulicSystem{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-25: add leaf hydraulic system
#     2022-May-27: move flow rates to a field FLOW
#     2022-Jun-13: use Union instead of Abstract... for type definition
#     2022-Jul-07: add e_crit as a field
#     2022-Jul-19: add dimension control to struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains leaf hydraulic system

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LeafHydraulics{FT<:AbstractFloat} <: AbstractHydraulicSystem{FT}
    # Dimensions
    "Dimension of xylem slices"
    DIM_XYLEM::Int = 5

    # General information of the hydraulic system
    "Leaf area `[m²]`"
    AREA::FT = 1500
    "Maximal extra-xylary hydraulic conductance `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    K_OX::FT = 100
    "Maximal leaf xylem hydraulic conductance per leaf area `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    K_SLA::FT = 0.04
    "Total capaciatance at Ψ = 0 `[mol m⁻²]`"
    V_MAXIMUM::FT = 20

    # Embedded structures
    "Flow profile"
    FLOW::Union{SteadyStateFlow{FT}, NonSteadyStateFlow{FT}} = SteadyStateFlow{FT}()
    "Pressure volume curve for storage"
    PVC::Union{LinearPVCurve{FT}, SegmentedPVCurve{FT}} = SegmentedPVCurve{FT}()
    "Vulnerability curve"
    VC::Union{LogisticVC{FT}, PowerVC{FT}, WeibullVC{FT}, ComplexVC{FT}} = WeibullVC{FT}()

    # Prognostic variables (used for ∂y∂t)
    "Vector of xylem water pressure history (normalized to 298.15 K) `[MPa]`"
    p_history::Vector{FT} = zeros(FT, DIM_XYLEM)
    "Current capaciatance at Ψ_leaf `[mol m⁻²]`"
    v_storage::FT = V_MAXIMUM

    # Prognostic variables (not used for ∂y∂t)
    "Leaf end water pressure at air-water interface `[MPa]`"
    p_leaf::FT = 0
    "Leaf xylem water pressure at the leaf base (upstream) `[MPa]`"
    p_ups::FT = 0

    # Cache variables
    "Critical flow rate `[mol s⁻¹ m⁻²]`"
    _e_crit::FT = 0
    "Vector of leaf kr history per element `[-]`"
    _k_history::Vector{FT} = ones(FT, DIM_XYLEM)
    "Leaf xylem water pressure at the downstream end of leaf xylem `[MPa]`"
    _p_dos::FT = 0
    "Vector of xylem water pressure `[MPa]`"
    _p_element::Vector{FT} = zeros(FT, DIM_XYLEM)
    "Pressure of storage"
    _p_storage::FT = 0
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-25: add root hydraulic system
#     2022-May-25: rename PV to PVC to be consistent with LeafHydraulics
#     2022-May-27: move flow rates to a field FLOW
#     2022-Jun-13: use Union instead of Abstract... for type definition
#     2022-Jul-19: add dimension control to struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains root hydraulic system

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct RootHydraulics{FT<:AbstractFloat} <: AbstractHydraulicSystem{FT}
    # Dimensions
    "Dimension of xylem slices"
    DIM_XYLEM::Int = 5

    # General information of the hydraulic system
    "Root cross-section area `[m²]`"
    AREA::FT = 1
    "Rhizosphere  conductance `[mol s⁻¹ MPa⁻¹]`"
    K_RHIZ::FT = 5e14
    "Maximal xylem hydraulic conductivity `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    K_X::FT = 25
    "Length `[m]`"
    L::FT = 1
    "Maximal storage per element `[mol]`"
    V_MAXIMUM::Vector{FT} = AREA * L / DIM_XYLEM * 6000 * ones(FT, DIM_XYLEM)
    "Root z difference `[m]`"
    ΔH::FT = 1

    # Embedded structures
    "Flow profile"
    FLOW::Union{SteadyStateFlow{FT}, NonSteadyStateFlow{FT}} = SteadyStateFlow{FT}()
    "Pressure volume curve for storage"
    PVC::Union{LinearPVCurve{FT}, SegmentedPVCurve{FT}} = LinearPVCurve{FT}()
    "Soil hydraulics"
    SH::Union{BrooksCorey{FT}, VanGenuchten{FT}} = VanGenuchten{FT}("Loam")
    "Vulnerability curve"
    VC::Union{LogisticVC{FT}, PowerVC{FT}, WeibullVC{FT}, ComplexVC{FT}} = WeibullVC{FT}()

    # Prognostic variables (used for ∂y∂t)
    "Vector of xylem water pressure history (normalized to 298.15 K) `[MPa]`"
    p_history::Vector{FT} = zeros(FT, DIM_XYLEM)
    "Storage per element `[mol]`"
    v_storage::Vector{FT} = V_MAXIMUM .* 1

    # Prognostic variables (not used for ∂y∂t)
    "Xylem water pressure at the downstream end of xylem `[MPa]`"
    p_dos::FT = 0
    "Soil matrix potential `[MPa]`"
    p_ups::FT = 0
    "Soil osmotic potential at 298.15 K `[MPa]"
    ψ_osm::FT = 0

    # Cache variables
    "Vector of leaf kr history per element"
    _k_history::Vector{FT} = ones(FT, DIM_XYLEM)
    "Vector of xylem water pressure `[MPa]`"
    _p_element::Vector{FT} = zeros(FT, DIM_XYLEM)
    "Xylem-rhizosphere interface water pressure `[MPa]`"
    _p_rhiz::FT = 0
    "Pressure of storage per element `[MPa]`"
    _p_storage::Vector{FT} = zeros(FT, DIM_XYLEM)
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-25: add stem hydraulic system
#     2022-May-25: rename PV to PVC to be consistent with LeafHydraulics
#     2022-May-27: move flow rates to a field FLOW
#     2022-Jun-13: use Union instead of Abstract... for type definition
#     2022-Jul-19: add dimension control to struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains stem hydraulic system

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct StemHydraulics{FT<:AbstractFloat} <: AbstractHydraulicSystem{FT}
    # Dimensions
    "Dimension of xylem slices"
    DIM_XYLEM::Int = 5

    # General information of the hydraulic system
    "Root cross-section area `[m²]`"
    AREA::FT = 1
    "Maximal xylem hydraulic conductivity (per root depth) `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    K_X::FT = 25
    "Length `[m]`"
    L::FT = 1
    "Maximal storage per element `[mol]`"
    V_MAXIMUM::Vector{FT} = AREA * L / DIM_XYLEM * 6000 * ones(FT, DIM_XYLEM)
    "Root z difference `[m]`"
    ΔH::FT = 1

    # Embedded structures
    "Flow profile"
    FLOW::Union{SteadyStateFlow{FT}, NonSteadyStateFlow{FT}} = SteadyStateFlow{FT}()
    "Pressure volume curve for storage"
    PVC::Union{LinearPVCurve{FT}, SegmentedPVCurve{FT}} = LinearPVCurve{FT}()
    "Vulnerability curve"
    VC::Union{LogisticVC{FT}, PowerVC{FT}, WeibullVC{FT}, ComplexVC{FT}} = WeibullVC{FT}()

    # Prognostic variables (used for ∂y∂t)
    "Vector of xylem water pressure history (normalized to 298.15 K) `[MPa]`"
    p_history::Vector{FT} = zeros(FT, DIM_XYLEM)
    "Storage per element `[mol]`"
    v_storage::Vector{FT} = V_MAXIMUM .* 1

    # Prognostic variables (not used for ∂y∂t)
    "Xylem water pressure at the downstream end of xylem `[MPa]`"
    p_dos::FT = 0
    "Soil matrix potential `[MPa]`"
    p_ups::FT = 0

    # Cache variables
    "Vector of leaf kr history per element"
    _k_history::Vector{FT} = ones(FT, DIM_XYLEM)
    "Vector of xylem water pressure `[MPa]`"
    _p_element::Vector{FT} = zeros(FT, DIM_XYLEM)
    "Pressure of storage per element"
    _p_storage::Vector{FT} = zeros(FT, DIM_XYLEM)
end
