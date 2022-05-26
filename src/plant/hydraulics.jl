#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-May-24: add abstract type for hydraulic organ
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierachy of AbstractHydraulicSystem:
"""
abstract type AbstractHydraulicSystem{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-May-25: add leaf hydraulic system
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
    "Maximal extra-xylary hydraulic conductance `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    K_OX::FT
    "Maximal leaf xylem hydraulic conductance per leaf area `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    K_SLA::FT
    "Number of xylem slices"
    N::Int
    "Pressure volume curve for storage"
    PVC::AbstractPVCurve{FT}
    "Total capaciatance at Ψ = 0 `[mol m⁻²]`"
    V_MAXIMUM::FT
    "Vulnerability curve"
    VC::AbstractXylemVC{FT}

    # prognostic variables that change with time
    "Flow rate in the xylem `[mol s⁻¹]`"
    flow::FT
    "Leaf xylem water pressure at the downstream end of leaf xylem `[MPa]`"
    p_dos::FT
    "Leaf end water pressure at air-water interface `[MPa]`"
    p_leaf::FT
    "Pressure of storage"
    p_storage::FT
    "Leaf xylem water pressure at the leaf base (upstream) `[MPa]`"
    p_ups::FT
    "Flow rate into the tissue (used for non-steady state) `[mol m⁻² s⁻¹]`"
    q_in::FT
    "Flow rate out of the tissue (used for non-steady state) `[mol m⁻² s⁻¹]`"
    q_out::FT
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
#     2022-May-25: add leaf hydraulic constructor
#
#######################################################################################################################################################################################################
"""

    LeafHydraulics{FT}(N::Int = 5; area::Number = 1500, k_ox::Number = 100, k_sla::Number = 0.04, v_max::Number = 20) where {FT<:AbstractFloat}

Constructor for leaf hydraulic system, given
- `N` Number of xylem slices in the system, default is 5
- `area` Leaf area
- `k_ox` Maximum extraxylary hydraulic conductance per leaf area
- `k_sla` Maximum leaf xylem hydraulic conductance per leaf area
- `v_max` Total water capacitance at Ψ = 0 per leaf area

---
# Examples
```julia
lhs = LeafHydraulics{Float64}();
lhs = LeafHydraulics{Float64}(N = 5);
lhs = LeafHydraulics{Float64}(N = 5; area = 20, k_ox = 50, k_sla = 0.1, v_max = 20);
```
"""
LeafHydraulics{FT}(N::Int = 5; area::Number = 1500, k_ox::Number = 100, k_sla::Number = 0.04, v_max::Number = 20) where {FT<:AbstractFloat} = (
    return LeafHydraulics{FT}(
                area,                   # AREA
                k_ox,                   # K_OX
                k_sla,                  # K_SLA
                N,                      # N
                SegmentedPVCurve{FT}(), # PVC
                v_max,                  # V_MAXIMUM
                WeibullVC{FT}(2,5),     # VC
                0,                      # flow
                0,                      # p_dos
                0,                      # p_leaf
                0,                      # p_storage
                0,                      # p_ups
                0,                      # q_in
                0,                      # q_out
                v_max,                  # v_storage
                ones(FT,N),             # k_history
                zeros(FT,N),            # p_element
                zeros(FT,N)             # p_history
    )
);
