###############################################################################
#
# Hydraulic system
#
###############################################################################
#= AbstractHydraulicSystem type tree
---> LeafHydraulics
---> RootHydraulics
---> StemHydraulics
=#
"""
    type AbstractHydraulicSystem

Hierarchy of AbstractHydraulicSystem
- [`LeafHydraulics`](@ref)
- [`StemHydraulics`](@ref)
"""
abstract type AbstractHydraulicSystem end




"""
    struct LeafHydraulics{FT<:AbstractFloat}

A struct that contains leaf hydraulics information.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct LeafHydraulics{FT<:AbstractFloat} <: AbstractHydraulicSystem
    # leaf hydraulic parameters
    "Weibull function (`k = k_max * exp( -(-p/B)^C )`) parameter B `[MPa]`"
    b    ::FT = FT(2.0)
    "Weibull function (`k = k_max * exp( -(-p/B)^C )`) parameter C"
    c    ::FT = FT(5.0)
    "Flow rate in the xylem `[mol s⁻¹]`"
    flow ::FT = FT(0)
    "Critical xylem pressure `[MPa]`"
    p_crt::FT = -b * log(FT(1e6)) ^ (1/c)
    "Maximal leaf hydraulic conductance per leaf area `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    k_sla::FT = FT(0.1)
    "Maximal extra-xylary hydraulic conductance `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    k_ox ::FT = FT(100)

    # flows and pressures (need to be updated with time)
    "Leaf xylem water pressure at the leaf base (upstream) `[MPa]`"
    p_ups ::FT = FT(0.0)
    "Leaf xylem water pressure at the downstream end of leaf xylem `[MPa]`"
    p_dos ::FT = FT(0.0)
    "Leaf end water pressure at air-water interface `[MPa]`"
    p_leaf::FT = FT(0.0)

    # pressure, k, and p_history profile
    "List of leaf k_max per element `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    k_element::Array{FT,1} =  ones(FT,10) .* 10 .* k_sla
    "List of leaf kr history per element"
    k_history::Array{FT,1} =  ones(FT,10)
    "List of xylem water pressure `[MPa]`"
    p_element::Array{FT,1} = zeros(FT,10)
    "List of xylem water pressure history (normalized to 298.15 K) `[MPa]`"
    p_history::Array{FT,1} = zeros(FT,10)

    # temperature related (uniform leaf temperature)
    "Relative surface tension"
    f_st ::FT = FT(1)
    "Relative viscosity"
    f_vis::FT = FT(1)
end




"""
    struct StemHydraulics{FT<:AbstractFloat}

A struct that contains stem hydraulics information.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct StemHydraulics{FT<:AbstractFloat} <: AbstractHydraulicSystem
    # stem hydraulic parameters
    "Weibull function (`k = k_max * exp( -(-p/B)^C )`) parameter B `[MPa]`"
    b    ::FT = FT(2.0)
    "Weibull function (`k = k_max * exp( -(-p/B)^C )`) parameter C"
    c    ::FT = FT(5.0)
    "Maximal leaf hydraulic conductance per leaf area `[mol s⁻¹ MPa⁻¹]`"
    k_max::FT = FT(0.1)
    "Flow rate in the xylem `[mol s⁻¹]`"
    flow ::FT = FT(0)

    # flows and pressures (need to be updated with time)
    "Xylem water pressure at the leaf base (upstream) `[MPa]`"
    p_ups ::FT = FT(0.0)
    "Xylem water pressure at the downstream end of leaf xylem `[MPa]`"
    p_dos ::FT = FT(0.0)

    # pressure, k, and p_history profile
    "List of leaf k_max per element `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    k_element::Array{FT,1} =  ones(FT,10) .* 10 .* k_max
    "List of leaf kr history per element"
    k_history::Array{FT,1} =  ones(FT,10)
    "List of xylem water pressure `[MPa]`"
    p_element::Array{FT,1} = zeros(FT,10)
    "List of xylem water pressure history (normalized to 298.15 K) `[MPa]`"
    p_history::Array{FT,1} = zeros(FT,10)

    # temperature related (uniform leaf temperature)
    "Relative surface tension"
    f_st ::FT = FT(1)
    "Relative viscosity"
    f_vis::FT = FT(1)
end
