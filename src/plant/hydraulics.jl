#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-May-24: add abstract type for hydraulic organ
#     2022-May-25: add documentation for type hierarchy
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractHydraulicSystem:
- [`LeafHydraulics`](@ref)
- [`RootHydraulics`](@ref)
"""
abstract type AbstractHydraulicSystem{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-May-25: add leaf hydraulic system
#     2022-May-27: move flow rates to a field FLOW
#     2022-Jun-13: use Union instead of Abstract... for type definition
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains leaf hydraulic system

# Fields

$(TYPEDFIELDS)

"""
mutable struct LeafHydraulics{FT} <: AbstractHydraulicSystem{FT}
    # parameters that do not change with time
    "Leaf area"
    AREA::FT
    "Flow profile"
    FLOW::Union{SteadyStateFlow{FT}, NonSteadyStateFlow{FT}}
    "Maximal extra-xylary hydraulic conductance `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    K_OX::FT
    "Maximal leaf xylem hydraulic conductance per leaf area `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    K_SLA::FT
    "Number of xylem slices"
    N::Int
    "Pressure volume curve for storage"
    PVC::Union{LinearPVCurve{FT}, SegmentedPVCurve{FT}}
    "Total capaciatance at Ψ = 0 `[mol m⁻²]`"
    V_MAXIMUM::FT
    "Vulnerability curve"
    VC::Union{LogisticVC{FT}, PowerVC{FT}, WeibullVC{FT}, ComplexVC{FT}}

    # prognostic variables that change with time
    "Leaf xylem water pressure at the downstream end of leaf xylem `[MPa]`"
    p_dos::FT
    "Leaf end water pressure at air-water interface `[MPa]`"
    p_leaf::FT
    "Pressure of storage"
    p_storage::FT
    "Leaf xylem water pressure at the leaf base (upstream) `[MPa]`"
    p_ups::FT
    "Current capaciatance at Ψ_leaf `[mol m⁻²]`"
    v_storage::FT

    # dignostic variables that change with time
    "Vector of leaf kr history per element"
    k_history::Vector{FT}
    "Vector of xylem water pressure `[MPa]`"
    p_element::Vector{FT}
    "Vector of xylem water pressure history (normalized to 298.15 K) `[MPa]`"
    p_history::Vector{FT}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-May-25: add leaf hydraulics constructor
#     2022-May-27: move flow rates to a field FLOW
#     2022-May-31: rename input variable name from steadystate to ssm
#
#######################################################################################################################################################################################################
"""

    LeafHydraulics{FT}(N::Int = 5; area::Number = 1500, k_ox::Number = 100, k_sla::Number = 0.04, v_max::Number = 20, ssm::Bool = true) where {FT<:AbstractFloat}

Constructor for leaf hydraulic system, given
- `N` Number of xylem slices in the system, default is 5
- `area` Leaf area
- `k_ox` Maximum extraxylary hydraulic conductance per leaf area
- `k_sla` Maximum leaf xylem hydraulic conductance per leaf area
- `v_max` Total water capacitance at Ψ = 0 per leaf area
- `ssm` Whether the flow rate is at steady state

---
# Examples
```julia
lhs = LeafHydraulics{Float64}();
lhs = LeafHydraulics{Float64}(N = 5);
lhs = LeafHydraulics{Float64}(N = 5; area = 20, k_ox = 50, k_sla = 0.1, v_max = 20, ssm = true);
```
"""
LeafHydraulics{FT}(N::Int = 5; area::Number = 1500, k_ox::Number = 100, k_sla::Number = 0.04, v_max::Number = 20, ssm::Bool = true) where {FT<:AbstractFloat} = (
    _flow = (ssm ? SteadyStateFlow{FT}(0) : NonSteadyStateFlow{FT}(N, true));

    return LeafHydraulics{FT}(
                area,                   # AREA
                _flow,                  # FLOW
                k_ox,                   # K_OX
                k_sla,                  # K_SLA
                N,                      # N
                SegmentedPVCurve{FT}(), # PVC
                v_max,                  # V_MAXIMUM
                WeibullVC{FT}(2,5),     # VC
                0,                      # p_dos
                0,                      # p_leaf
                0,                      # p_storage
                0,                      # p_ups
                v_max,                  # v_storage
                ones(FT,N),             # k_history
                zeros(FT,N),            # p_element
                zeros(FT,N)             # p_history
    )
);


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-May-25: add root hydraulic system
#     2022-May-25: rename PV to PVC to be consistent with LeafHydraulics
#     2022-May-27: move flow rates to a field FLOW
#     2022-Jun-13: use Union instead of Abstract... for type definition
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains root hydraulic system

# Fields

$(TYPEDFIELDS)

"""
mutable struct RootHydraulics{FT} <: AbstractHydraulicSystem{FT}
    # parameters that do not change with time
    "Root cross-section area `[m²]`"
    AREA::FT
    "Flow profile"
    FLOW::Union{SteadyStateFlow{FT}, NonSteadyStateFlow{FT}}
    "Maximal hydraulic conductance `[mol s⁻¹ MPa⁻¹]`"
    K_MAX::FT
    "Rhizosphere  conductance `[mol s⁻¹ MPa⁻¹]`"
    K_RHIZ::FT
    "Maximal xylem hydraulic conductivity `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    K_X::FT
    "Length `[m]`"
    L::FT
    "Number of xylem slices"
    N::Int
    "Pressure volume curve for storage"
    PVC::Union{LinearPVCurve{FT}, SegmentedPVCurve{FT}}
    "Soil hydraulics"
    SH::Union{BrooksCorey{FT}, VanGenuchten{FT}}
    "Maximal storage per element `[mol]`"
    V_MAXIMUM::Vector{FT}
    "Vulnerability curve"
    VC::Union{LogisticVC{FT}, PowerVC{FT}, WeibullVC{FT}, ComplexVC{FT}}
    "Root z difference `[m]`"
    ΔH::FT

    # prognostic variables that change with time
    "Xylem water pressure at the downstream end of xylem `[MPa]`"
    p_dos::FT
    "Xylem-rhizosphere interface water pressure `[MPa]`"
    p_rhiz::FT
    "Pressure of storage per element"
    p_storage::Vector{FT}
    "Soil matrix potential `[MPa]`"
    p_ups::FT
    "Storage per element `[mol]`"
    v_storage::Vector{FT}
    "Soil osmotic potential at 298.15 K `[MPa]"
    ψ_osm::FT

    # dignostic variables that change with time
    "Vector of leaf kr history per element"
    k_history::Vector{FT}
    "Vector of xylem water pressure `[MPa]`"
    p_element::Vector{FT}
    "Vector of xylem water pressure history (normalized to 298.15 K) `[MPa]`"
    p_history::Vector{FT}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-May-25: add root hydraulics constructor
#     2022-May-27: move flow rates to a field FLOW
#     2022-May-31: rename input variable name from steadystate to ssm
# Bug fixes:
#     2022-May-31: set the isleaf option to false
#
#######################################################################################################################################################################################################
"""

    RootHydraulics{FT}(N::Int = 5; area::Number = 1, k_x::Number = 25, Δh::Number = 1, Δl::Number = 1, ssm::Bool = true) where {FT<:AbstractFloat}

Constructor for root hydraulic system, given
- `N` Number of xylem slices in the system, default is 5
- `area` Root crosssection area
- `k_x` Maximum root xylem hydraulic conductivity per crosssection area per root length
- `Δh` Root depth
- `Δl` Root length
- `ssm` Whether the flow rate is at steady state

---
# Examples
```julia
rhs = RootHydraulics{Float64}();
rhs = RootHydraulics{Float64}(N = 5);
rhs = RootHydraulics{Float64}(N = 5; area = 1, k_x = 50, Δh = 1, Δl = 2, ssm = true);
```
"""
RootHydraulics{FT}(N::Int = 5; area::Number = 1, k_x::Number = 25, Δh::Number = 1, Δl::Number = 1, ssm::Bool = true) where {FT<:AbstractFloat} = (
    _flow = (ssm ? SteadyStateFlow{FT}(0) : NonSteadyStateFlow{FT}(N, false));

    return RootHydraulics{FT}(
        area,                               # AREA
        _flow,                              # FLOW
        k_x * area / Δl,                    # K_MAX
        5e14,                               # K_RHIZ
        k_x,                                # K_X
        Δl,                                 # L
        N,                                  # N
        LinearPVCurve{FT}(),                # PVC
        VanGenuchten{FT}("Silt"),           # SH
        area * Δh / N * 6000 * ones(FT,N),  # V_MAXIMUM
        WeibullVC{FT}(2,5),                 # VC
        Δh,                                 # ΔH
        0,                                  # p_dos
        0,                                  # p_rhiz
        zeros(FT,N),                        # p_storage
        0,                                  # p_ups
        zeros(FT,N),                        # v_storage
        0,                                  # ψ_osm
        ones(FT,N),                         # k_history
        zeros(FT,N),                        # p_element
        zeros(FT,N)                         # p_history
    )
);


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-May-25: add stem hydraulic system
#     2022-May-25: rename PV to PVC to be consistent with LeafHydraulics
#     2022-May-27: move flow rates to a field FLOW
#     2022-Jun-13: use Union instead of Abstract... for type definition
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains stem hydraulic system

# Fields

$(TYPEDFIELDS)

"""
mutable struct StemHydraulics{FT} <: AbstractHydraulicSystem{FT}
    # parameters that do not change with time
    "Root cross-section area `[m²]`"
    AREA::FT
    "Flow profile"
    FLOW::Union{SteadyStateFlow{FT}, NonSteadyStateFlow{FT}}
    "Maximal hydraulic conductance `[mol s⁻¹ MPa⁻¹]`"
    K_MAX::FT
    "Maximal xylem hydraulic conductivity (per root depth) `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    K_X::FT
    "Length `[m]`"
    L::FT
    "Number of xylem slices"
    N::Int
    "Pressure volume curve for storage"
    PVC::Union{LinearPVCurve{FT}, SegmentedPVCurve{FT}}
    "Maximal storage per element `[mol]`"
    V_MAXIMUM::Vector{FT}
    "Vulnerability curve"
    VC::Union{LogisticVC{FT}, PowerVC{FT}, WeibullVC{FT}, ComplexVC{FT}}
    "Root z difference `[m]`"
    ΔH::FT

    # prognostic variables that change with time
    "Xylem water pressure at the downstream end of xylem `[MPa]`"
    p_dos::FT
    "Pressure of storage per element"
    p_storage::Vector{FT}
    "Soil matrix potential `[MPa]`"
    p_ups::FT
    "Storage per element `[mol]`"
    v_storage::Vector{FT}

    # dignostic variables that change with time
    "Vector of leaf kr history per element"
    k_history::Vector{FT}
    "Vector of xylem water pressure `[MPa]`"
    p_element::Vector{FT}
    "Vector of xylem water pressure history (normalized to 298.15 K) `[MPa]`"
    p_history::Vector{FT}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-May-25: add stem hydraulics constructor
#     2022-May-27: move flow rates to a field FLOW
#     2022-May-31: rename input variable name from steadystate to ssm
# Bug fixes:
#     2022-May-31: set the isleaf option to false
#
#######################################################################################################################################################################################################
"""

    StemHydraulics{FT}(N::Int = 5; area::Number = 1, k_x::Number = 25, Δh::Number = 1, Δl::Number = 1, ssm::Bool = true) where {FT<:AbstractFloat}

Constructor for stem hydraulic system, given
- `N` Number of xylem slices in the system, default is 5
- `area` Root crosssection area
- `k_x` Maximum stem xylem hydraulic conductivity per crosssection area per stem length
- `Δh` Stem height
- `Δl` Stem length
- `ssm` Whether the flow rate is at steady state

---
# Examples
```julia
rhs = StemHydraulics{Float64}();
rhs = StemHydraulics{Float64}(N = 5);
rhs = StemHydraulics{Float64}(N = 5; area = 1, k_x = 50, Δh = 2, ssm = true);
```
"""
StemHydraulics{FT}(N::Int = 5; area::Number = 1, k_x::Number = 25, Δh::Number = 1, Δl::Number = 1, ssm::Bool = true) where {FT<:AbstractFloat} = (
    _flow = (ssm ? SteadyStateFlow{FT}(0) : NonSteadyStateFlow{FT}(N, false));

    return StemHydraulics{FT}(
        area,                               # AREA
        _flow,                              # FLOW
        k_x * area / Δl,                    # K_MAX
        k_x,                                # K_X
        Δl,                                 # L
        N,                                  # N
        LinearPVCurve{FT}(),                # PVC
        area * Δh / N * 6000 * ones(FT,N),  # V_MAXIMUM
        WeibullVC{FT}(2,5),                 # VC
        Δh,                                 # ΔH
        0,                                  # p_dos
        zeros(FT,N),                        # p_storage
        0,                                  # p_ups
        zeros(FT,N),                        # v_storage
        ones(FT,N),                         # k_history
        zeros(FT,N),                        # p_element
        zeros(FT,N)                         # p_history
    )
);
