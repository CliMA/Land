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
# Changes to this type
# General
#     2022-May-25: add leaf hydraulic system
#     2022-May-27: move flow rates to a field FLOW
#     2022-Jun-13: use Union instead of Abstract... for type definition
#     2022-Jul-07: add e_crit as a field
#     2022-Jul-18: use kwdef for the constructor
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
    # dimensions
    "Dimension of xylem slices"
    DIM_XYLEM::Int = 5

    # parameters that do not change with time
    "Leaf area"
    AREA::FT = 1500
    "Flow profile"
    FLOW::Union{SteadyStateFlow{FT}, NonSteadyStateFlow{FT}} = SteadyStateFlow{FT}()
    "Maximal extra-xylary hydraulic conductance `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    K_OX::FT = 100
    "Maximal leaf xylem hydraulic conductance per leaf area `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    K_SLA::FT = 0.04
    "Pressure volume curve for storage"
    PVC::Union{LinearPVCurve{FT}, SegmentedPVCurve{FT}} = SegmentedPVCurve{FT}()
    "Total capaciatance at Ψ = 0 `[mol m⁻²]`"
    V_MAXIMUM::FT = 20
    "Vulnerability curve"
    VC::Union{LogisticVC{FT}, PowerVC{FT}, WeibullVC{FT}, ComplexVC{FT}} = WeibullVC{FT}()

    # prognostic variables that change with time
    "Leaf xylem water pressure at the downstream end of leaf xylem `[MPa]`"
    p_dos::FT = 0
    "Leaf end water pressure at air-water interface `[MPa]`"
    p_leaf::FT = 0
    "Pressure of storage"
    p_storage::FT = 0
    "Leaf xylem water pressure at the leaf base (upstream) `[MPa]`"
    p_ups::FT = 0
    "Current capaciatance at Ψ_leaf `[mol m⁻²]`"
    v_storage::FT = V_MAXIMUM

    # dignostic variables that change with time
    "Critical flow rate `[mol s⁻¹ m⁻²]`"
    e_crit::FT = 0
    "Vector of leaf kr history per element"
    k_history::Vector{FT} = ones(FT, DIM_XYLEM)
    "Vector of xylem water pressure `[MPa]`"
    p_element::Vector{FT} = zeros(FT, DIM_XYLEM)
    "Vector of xylem water pressure history (normalized to 298.15 K) `[MPa]`"
    p_history::Vector{FT} = zeros(FT, DIM_XYLEM)
end


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-May-25: add root hydraulic system
#     2022-May-25: rename PV to PVC to be consistent with LeafHydraulics
#     2022-May-27: move flow rates to a field FLOW
#     2022-Jun-13: use Union instead of Abstract... for type definition
#     2022-Jul-18: use kwdef for the constructor
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
    # dimensions
    "Dimension of xylem slices"
    DIM_XYLEM::Int = 5

    # parameters that do not change with time
    "Root cross-section area `[m²]`"
    AREA::FT = 1
    "Flow profile"
    FLOW::Union{SteadyStateFlow{FT}, NonSteadyStateFlow{FT}} = SteadyStateFlow{FT}()
    "Rhizosphere  conductance `[mol s⁻¹ MPa⁻¹]`"
    K_RHIZ::FT = 5e14
    "Maximal xylem hydraulic conductivity `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    K_X::FT = 25
    "Length `[m]`"
    L::FT = 1
    "Pressure volume curve for storage"
    PVC::Union{LinearPVCurve{FT}, SegmentedPVCurve{FT}} = LinearPVCurve{FT}()
    "Soil hydraulics"
    SH::Union{BrooksCorey{FT}, VanGenuchten{FT}} = VanGenuchten{FT}("Loam")
    "Maximal storage per element `[mol]`"
    V_MAXIMUM::Vector{FT} = AREA * L / DIM_XYLEM * 6000 * ones(FT, DIM_XYLEM)
    "Vulnerability curve"
    VC::Union{LogisticVC{FT}, PowerVC{FT}, WeibullVC{FT}, ComplexVC{FT}} = WeibullVC{FT}()
    "Root z difference `[m]`"
    ΔH::FT = 1

    # dignostic variables that change with time
    "Vector of leaf kr history per element"
    k_history::Vector{FT} = ones(FT, DIM_XYLEM)
    "Xylem water pressure at the downstream end of xylem `[MPa]`"
    p_dos::FT = 0
    "Vector of xylem water pressure `[MPa]`"
    p_element::Vector{FT} = zeros(FT, DIM_XYLEM)
    "Vector of xylem water pressure history (normalized to 298.15 K) `[MPa]`"
    p_history::Vector{FT} = zeros(FT, DIM_XYLEM)
    "Xylem-rhizosphere interface water pressure `[MPa]`"
    p_rhiz::FT = 0
    "Pressure of storage per element"
    p_storage::Vector{FT} = zeros(FT, DIM_XYLEM)
    "Soil matrix potential `[MPa]`"
    p_ups::FT = 0
    "Storage per element `[mol]`"
    v_storage::Vector{FT} = V_MAXIMUM .* 1
    "Soil osmotic potential at 298.15 K `[MPa]"
    ψ_osm::FT = 0
end


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-May-25: add stem hydraulic system
#     2022-May-25: rename PV to PVC to be consistent with LeafHydraulics
#     2022-May-27: move flow rates to a field FLOW
#     2022-Jun-13: use Union instead of Abstract... for type definition
#     2022-Jul-18: use kwdef for the constructor
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
    # dimensions
    "Dimension of xylem slices"
    DIM_XYLEM::Int = 5

    # parameters that do not change with time
    "Root cross-section area `[m²]`"
    AREA::FT = 1
    "Flow profile"
    FLOW::Union{SteadyStateFlow{FT}, NonSteadyStateFlow{FT}} = SteadyStateFlow{FT}()
    "Maximal xylem hydraulic conductivity (per root depth) `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    K_X::FT = 25
    "Length `[m]`"
    L::FT = 1
    "Pressure volume curve for storage"
    PVC::Union{LinearPVCurve{FT}, SegmentedPVCurve{FT}} = LinearPVCurve{FT}()
    "Maximal storage per element `[mol]`"
    V_MAXIMUM::Vector{FT} = AREA * L / DIM_XYLEM * 6000 * ones(FT, DIM_XYLEM)
    "Vulnerability curve"
    VC::Union{LogisticVC{FT}, PowerVC{FT}, WeibullVC{FT}, ComplexVC{FT}} = WeibullVC{FT}()
    "Root z difference `[m]`"
    ΔH::FT = 1

    # dignostic variables that change with time
    "Vector of leaf kr history per element"
    k_history::Vector{FT} = ones(FT, DIM_XYLEM)
    "Xylem water pressure at the downstream end of xylem `[MPa]`"
    p_dos::FT = 0
    "Vector of xylem water pressure `[MPa]`"
    p_element::Vector{FT} = zeros(FT, DIM_XYLEM)
    "Vector of xylem water pressure history (normalized to 298.15 K) `[MPa]`"
    p_history::Vector{FT} = zeros(FT, DIM_XYLEM)
    "Pressure of storage per element"
    p_storage::Vector{FT} = zeros(FT, DIM_XYLEM)
    "Soil matrix potential `[MPa]`"
    p_ups::FT = 0
    "Storage per element `[mol]`"
    v_storage::Vector{FT} = V_MAXIMUM .* 1
end
