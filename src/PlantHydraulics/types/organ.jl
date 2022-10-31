###############################################################################
#
# Segmented Hydraulic system
#
###############################################################################
#= AbstractHydraulicOrgan type tree
---> LeafHydraulics
---> RootHydraulics
---> StemHydraulics
=#
"""
    abstract type AbstractHydraulicOrgan{FT}

Hierarchy of AbstractHydraulicOrgan
- [`LeafHydraulics`](@ref)
- [`RootHydraulics`](@ref)
- [`StemHydraulics`](@ref)
"""
abstract type AbstractHydraulicOrgan{FT<:AbstractFloat} end




"""
    mutable struct LeafHydraulics{FT}

A struct that contains leaf hydraulics information.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct LeafHydraulics{FT} <: AbstractHydraulicOrgan{FT}
    # number of xylem slices
    N::Int = 5

    # leaf hydraulic parameters
    "Leaf area `[m²]`"
    area ::FT = FT(1500)
    "Maximal extra-xylary hydraulic conductance `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    k_ox ::FT = FT(100)
    "Maximal leaf hydraulic conductance per leaf area `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    k_sla::FT = FT(0.04)
    "Vulnerability curve"
    vc::AbstractXylemVC{FT} = WeibullSingle{FT}()
    "Critical xylem pressure `[MPa]`"
    p_crt::FT = -vc.b * log(FT(1000)) ^ (1/vc.c)

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
    k_element::Array{FT,1} =  ones(FT,N) .* N .* k_sla
    "List of leaf kr history per element"
    k_history::Array{FT,1} =  ones(FT,N)
    "List of xylem water pressure `[MPa]`"
    p_element::Array{FT,1} = zeros(FT,N)
    "List of xylem water pressure history (normalized to 298.15 K) `[MPa]`"
    p_history::Array{FT,1} = zeros(FT,N)

    # temperature related (uniform leaf temperature), update with time
    "Relative surface tension"
    f_st ::FT = FT(1)
    "Relative viscosity"
    f_vis::FT = FT(1)
    "Temperature memory `[K]`"
    T_old::FT = T₂₅(FT)
    "Upstream sap temperature `[K]`"
    T_sap::FT = T₂₅(FT)

    # capacitance
    "Pressure volume curve for storage"
    pv::AbstractCapacity{FT} = PVCurveSegmented{FT}()
    "Pressure of storage"
    p_storage::FT = 0
    "Total capaciatance at Ψ = 0 `[mol m⁻²]`"
    v_maximum::FT = 20
    "Current capaciatance at Ψ_leaf `[mol m⁻²]`"
    v_storage::FT = 20
    "Flow rate into the tissue (used for non-steady state) `[mol m⁻² s⁻¹]`"
    q_in     ::FT = 0
    "Flow rate out of the tissue (used for non-steady state) `[mol m⁻² s⁻¹]`"
    q_out    ::FT = 0
end




"""
    mutable struct RootHydraulics{FT}

A struct that contains root hydraulics information.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct RootHydraulics{FT} <: AbstractHydraulicOrgan{FT}
    # number of xylem slices
    N::Int = 5

    # root hydraulic parameters
    "Root cross-section area `[m²]`"
    area ::FT = FT(1)
    "Maximal hydraulic conductance `[mol s⁻¹ MPa⁻¹]`"
    k_max::FT = FT(25)
    "Maximal xylem hydraulic conductivity `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    k_s  ::FT = FT(250)
    "Vulnerability curve"
    vc::AbstractXylemVC{FT} = WeibullSingle{FT}()
    "Root z difference `[m]`"
    Δh   ::FT = FT(1.0)

    # soil parameters
    "Rhizosphere  conductance `[mol s⁻¹ MPa⁻¹]`"
    k_rhiz   ::FT     = FT(5e14)
    "Soil hydraulics"
    sh::Union{BrooksCorey{FT}, VanGenuchten{FT}} = BrooksCorey{FT}()

    # flows and pressures (need to be updated with time)
    "Flow rate in the xylem `[mol s⁻¹]`"
    flow  ::FT = FT(0.0)
    "Xylem water pressure at the downstream end of xylem `[MPa]`"
    p_dos ::FT = FT(0.0)
    "Xylem-rhizosphere interface water pressure `[MPa]`"
    p_rhiz::FT = FT(0.0)
    "Soil matrix potential `[MPa]`"
    p_ups ::FT = FT(0.0)
    "Soil osmotic potential at 298.15 K `[MPa]"
    p_osm ::FT = FT(0.0)

    # pressure, k, and p_history profile
    "List of k_max per element `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    k_element::Array{FT,1} =  ones(FT,N) .* N .* k_max
    "List of kr history per element"
    k_history::Array{FT,1} =  ones(FT,N)
    "List of xylem water pressure `[MPa]`"
    p_element::Array{FT,1} = zeros(FT,N)
    "List of pressure drop caused by gravity `[MPa]`"
    p_gravity::Array{FT,1} = zeros(FT,N) .+ ρg_MPa(FT) .* Δh ./ N;
    "List of xylem water pressure history (normalized to 298.15 K) `[MPa]`"
    p_history::Array{FT,1} = zeros(FT,N)

    # temperature related (uniform temperature), update with time
    "Relative surface tension"
    f_st ::FT = FT(1)
    "Relative viscosity"
    f_vis::FT = FT(1)
    "Temperature memory `[K]`"
    T_old::FT = T₂₅(FT)
    "Upstream sap temperature `[K]`"
    T_sap::FT = T₂₅(FT)

    # capacitance
    "Pressure volume curve for storage"
    pv::AbstractCapacity{FT} = PVCurveLinear{FT}()
    "Pressure of storage per element"
    p_storage::Array{FT,1} = zeros(FT,N)
    "Maximal storage per element `[mol]`"
    v_maximum::Array{FT,1} = area * Δh / N * 6000 * ones(FT,N)
    "Storage per element `[mol]`"
    v_storage::Array{FT,1} = area * Δh / N * 6000 * ones(FT,N)
    "List of xylem water flow `[mol m⁻²]`"
    q_element::Array{FT,1} = zeros(FT,N)
    "List of buffer water flow `[mol m⁻²]`"
    q_buffer ::Array{FT,1} = zeros(FT,N)
    "List of diiferntial water flow `[mol m⁻²]`"
    q_diff   ::Array{FT,1} = zeros(FT,N)
    "Flow rate into the tissue (used for non-steady state) `[mol s⁻¹]`"
    q_in ::FT = 0
    "Flow rate out of the tissue (used for non-steady state) `[mol s⁻¹]`"
    q_out::FT = 0
end




"""
    mutable struct StemHydraulics{FT}

A struct that contains stem hydraulics information.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct StemHydraulics{FT} <: AbstractHydraulicOrgan{FT}
    # number of xylem slices
    N::Int = 5

    # stem hydraulic parameters
    "Stem cross-section area `[m²]`"
    area ::FT = FT(1)
    "Maximal hydraulic conductance `[mol s⁻¹ MPa⁻¹]`"
    k_max::FT = FT(50)
    "Maximal xylem hydraulic conductivity `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    k_s  ::FT = FT(250)
    "Vulnerability curve"
    vc::AbstractXylemVC{FT} = WeibullSingle{FT}()
    "Stem height difference `[m]`"
    Δh   ::FT = FT(5.0)

    # flows and pressures (need to be updated with time)
    "Flow rate in the xylem `[mol s⁻¹]`"
    flow  ::FT = FT(0.0)
    "Xylem water pressure at the downstream end of xylem `[MPa]`"
    p_dos::FT = FT(0.0)
    "Xylem water pressure at the base (upstream) `[MPa]`"
    p_ups::FT = FT(0.0)

    # pressure, k, and p_history profile
    "List of k_max per element `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    k_element::Array{FT,1} =  ones(FT,N) .* N .* k_max
    "List of kr history per element"
    k_history::Array{FT,1} =  ones(FT,N)
    "List of xylem water pressure `[MPa]`"
    p_element::Array{FT,1} = zeros(FT,N)
    "List of pressure drop caused by gravity `[MPa]`"
    p_gravity::Array{FT,1} = zeros(FT,N) .+ ρg_MPa(FT) .* Δh ./ N;
    "List of xylem water pressure history (normalized to 298.15 K) `[MPa]`"
    p_history::Array{FT,1} = zeros(FT,N)

    # temperature related (uniform temperature), update with time
    "Relative surface tension"
    f_st ::FT = FT(1)
    "Relative viscosity"
    f_vis::FT = FT(1)
    "Temperature memory `[K]`"
    T_old::FT = T₂₅(FT)
    "Upstream sap temperature `[K]`"
    T_sap::FT = T₂₅(FT)

    # capacitance
    "Pressure volume curve for storage"
    pv::AbstractCapacity{FT} = PVCurveLinear{FT}()
    "Pressure of storage per element"
    p_storage::Array{FT,1} = zeros(FT,N)
    "Maximal storage per element `[mol]`"
    v_maximum::Array{FT,1} = area * Δh / N * 6000 * ones(FT,N)
    "Storage per element `[mol]`"
    v_storage::Array{FT,1} = area * Δh / N * 6000 * ones(FT,N)
    "List of xylem water flow `[mol m⁻²]`"
    q_element::Array{FT,1} = zeros(FT,N)
    "Flow rate into the tissue (used for non-steady state) `[mol s⁻¹]`"
    q_in ::FT = 0
    "Flow rate out of the tissue (used for non-steady state) `[mol s⁻¹]`"
    q_out::FT = 0
end
