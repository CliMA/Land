###############################################################################
#
# Segmented Hydraulic system
#
###############################################################################
#= AbstractHydraulicSystem type tree
---> LeafHydraulics
---> RootHydraulics
---> StemHydraulics
=#
"""
    abstract type AbstractHydraulicSystem{FT}

Hierarchy of AbstractHydraulicSystem
- [`LeafHydraulics`](@ref)
- [`RootHydraulics`](@ref)
- [`StemHydraulics`](@ref)
"""
abstract type AbstractHydraulicSystem{FT} end




"""
    mutable struct LeafHydraulics{FT<:AbstractFloat}

A struct that contains leaf hydraulics information.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct LeafHydraulics{FT<:AbstractFloat} <: AbstractHydraulicSystem{FT}
    # leaf hydraulic parameters
    "Leaf area `[m²]`"
    area ::FT = FT(150)
    "Maximal extra-xylary hydraulic conductance `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    k_ox ::FT = FT(100)
    "Maximal leaf hydraulic conductance per leaf area `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    k_sla::FT = FT(0.1)
    "Vulnerability curve"
    vc::AbstractVulnerability{FT} = WeibullSingle{FT}()
    "Critical xylem pressure `[MPa]`"
    p_crt::FT = -vc.b * log(FT(1e6)) ^ (1/vc.c)

    # flows and pressures (need to be updated with time)
    "Flow rate in the xylem `[mol s⁻¹]`"
    flow  ::FT = FT(0)
    "Leaf xylem water pressure at the downstream end of leaf xylem `[MPa]`"
    p_dos ::FT = FT(0.0)
    "Leaf end water pressure at air-water interface `[MPa]`"
    p_leaf::FT = FT(0.0)
    "Leaf xylem water pressure at the leaf base (upstream) `[MPa]`"
    p_ups ::FT = FT(0.0)

    # pressure, k, and p_history profile
    "List of leaf k_max per element `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    k_element::Array{FT,1} =  ones(FT,10) .* 10 .* k_sla
    "List of leaf kr history per element"
    k_history::Array{FT,1} =  ones(FT,10)
    "List of xylem water pressure `[MPa]`"
    p_element::Array{FT,1} = zeros(FT,10)
    "List of xylem water pressure history (normalized to 298.15 K) `[MPa]`"
    p_history::Array{FT,1} = zeros(FT,10)

    # temperature related (uniform leaf temperature), update with time
    "Relative surface tension"
    f_st ::FT = FT(1)
    "Relative viscosity"
    f_vis::FT = FT(1)
    "Upstream sap temperature `[K]`"
    T_sap::FT = FT(K_25)
end




"""
    mutable struct RootHydraulics{FT<:AbstractFloat}

A struct that contains root hydraulics information.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct RootHydraulics{FT<:AbstractFloat} <: AbstractHydraulicSystem{FT}
    # root hydraulic parameters
    "Root cross-section area `[m²]`"
    area ::FT = FT(0.02)
    "Maximal hydraulic conductance `[mol s⁻¹ MPa⁻¹]`"
    k_max::FT = FT(0.5)
    "Maximal xylem hydraulic conductivity `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    k_s  ::FT = FT(250)
    "Vulnerability curve"
    vc::AbstractVulnerability{FT} = WeibullSingle{FT}()
    "Root z difference `[m]`"
    Δh   ::FT = FT(1.0)

    # soil parameters
    "Rhizosphere  conductance `[mol s⁻¹ MPa⁻¹]`"
    k_rhiz   ::FT     = FT(1e13)
    "Soil type"
    soil_type::String = "Sandy Clay Loam"
    "α is related to the inverse of the air entry suction, α > 0"
    soil_α   ::FT     = FT(602.0419)
    "Measure of the pore-size distribution"
    soil_n   ::FT     = FT(1.48)
    "1 - 1/n"
    soil_m   ::FT     = 1 - 1/soil_n

    # flows and pressures (need to be updated with time)
    "Flow rate in the xylem `[mol s⁻¹]`"
    flow  ::FT = FT(0.0)
    "Xylem water pressure at the downstream end of xylem `[MPa]`"
    p_dos ::FT = FT(0.0)
    "Xylem-rhizosphere interface water pressure `[MPa]`"
    p_rhiz::FT = FT(0.0)
    "Soil matrix potential `[MPa]`"
    p_ups ::FT = FT(0.0)

    # pressure, k, and p_history profile
    "List of k_max per element `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    k_element::Array{FT,1} =  ones(FT,10) .* 10 .* k_max
    "List of kr history per element"
    k_history::Array{FT,1} =  ones(FT,10)
    "List of xylem water pressure `[MPa]`"
    p_element::Array{FT,1} = zeros(FT,10)
    "List of pressure drop caused by gravity `[MPa]`"
    p_gravity::Array{FT,1} = zeros(FT,10) .+ ρg_MPa .* Δh ./ 10;
    "List of xylem water pressure history (normalized to 298.15 K) `[MPa]`"
    p_history::Array{FT,1} = zeros(FT,10)

    # temperature related (uniform temperature), update with time
    "Relative surface tension"
    f_st ::FT = FT(1)
    "Relative viscosity"
    f_vis::FT = FT(1)
    "Upstream sap temperature `[K]`"
    T_sap::FT = FT(K_25)
end




"""
    mutable struct StemHydraulics{FT<:AbstractFloat}

A struct that contains stem hydraulics information.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct StemHydraulics{FT<:AbstractFloat} <: AbstractHydraulicSystem{FT}
    # stem hydraulic parameters
    "Stem cross-section area `[m²]`"
    area ::FT = FT(0.1)
    "Maximal hydraulic conductance `[mol s⁻¹ MPa⁻¹]`"
    k_max::FT = FT(5.0)
    "Maximal xylem hydraulic conductivity `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    k_s  ::FT = FT(250)
    "Vulnerability curve"
    vc::AbstractVulnerability{FT} = WeibullSingle{FT}()
    "Stem height difference `[m]`"
    Δh   ::FT = FT(5.0)

    # flows and pressures (need to be updated with time)
    "Xylem water pressure at the downstream end of xylem `[MPa]`"
    p_dos::FT = FT(0.0)
    "Xylem water pressure at the base (upstream) `[MPa]`"
    p_ups::FT = FT(0.0)

    # pressure, k, and p_history profile
    "List of k_max per element `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    k_element::Array{FT,1} =  ones(FT,10) .* 10 .* k_max
    "List of kr history per element"
    k_history::Array{FT,1} =  ones(FT,10)
    "List of xylem water pressure `[MPa]`"
    p_element::Array{FT,1} = zeros(FT,10)
    "List of pressure drop caused by gravity `[MPa]`"
    p_gravity::Array{FT,1} = zeros(FT,10) .+ ρg_MPa .* Δh ./ 10;
    "List of xylem water pressure history (normalized to 298.15 K) `[MPa]`"
    p_history::Array{FT,1} = zeros(FT,10)

    # temperature related (uniform temperature), update with time
    "Relative surface tension"
    f_st ::FT = FT(1)
    "Relative viscosity"
    f_vis::FT = FT(1)
    "Upstream sap temperature `[K]`"
    T_sap::FT = FT(K_25)
end
