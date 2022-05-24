###############################################################################
#
# Capacitance for plants
#
###############################################################################
"""
    $(TYPEDEF)

Hierachy of AbstractCapacity
- [`PVCurveLinear`](@ref)
- [`PVCurveSegmented`](@ref)
"""
abstract type AbstractCapacity{FT<:AbstractFloat} end




"""
    $(TYPEDEF)

Struct that contains information for linear PV curve

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct PVCurveLinear{FT<:AbstractFloat} <: AbstractCapacity{FT}
    "Slope of the linear PV curve (relative to maximum) `[MPa⁻¹]`"
    slope   ::FT = 0.2
    "Conductance for refilling (relative to maximum) `[MPa⁻¹ s⁻¹]`"
    k_refill::FT = 0.0001
end




"""
    $(TYPEDEF)


"""
Base.@kwdef mutable struct PVCurveSegmented{FT} <: AbstractCapacity{FT}
    "n_o / maximum V `[mol m⁻³]`"
    c_all   ::FT = 300.0
    "Bulk modulus of elasticity `[MPa]`"
    ϵ_bulk  ::FT = 20.0
    "Apoplastic water content relative to maximum water volume"
    RWC_apo ::FT = 0.2
    "Relative water content at turgor loss point"
    RWC_TLP ::FT = 0.8
    "Conductance for refilling (relative to maximum) `[MPa⁻¹ s⁻¹]`"
    k_refill::FT = 0.0001
end
