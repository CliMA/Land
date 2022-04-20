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






#=

###############################################################################
#
# Vulnerability curves for xylem
#
###############################################################################
"""
    $(TYPEDEF)

Hierachy of AbstractXylemVC
- [`LogisticSingle`](@ref)
- [`PowerSingle`](@ref)
- [`WeibullSingle`](@ref)
- [`WeibullDual`](@ref)
"""
abstract type AbstractXylemVC{FT<:AbstractFloat} end




"""
    $(TYPEDEF)

Struct that contains single modified logistc function parameters.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct LogisticSingle{FT} <: AbstractXylemVC{FT}
    "Multiplier to exp"
    a::FT = 2
    "Multiplier to P `[MPa⁻¹]`"
    b::FT = 2
end




"""
    $(TYPEDEF)

Struct that contains single power function parameters.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct PowerSingle{FT} <: AbstractXylemVC{FT}
    "Multiplier to power function, `[MPa⁻ᵇ]`"
    a::FT = 2
    "Power to pressure"
    b::FT = 2
end




"""
    $(TYPEDEF)

Struct that contains dual Weibull function parameters.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct WeibullDual{FT} <: AbstractXylemVC{FT}
    "B of first part `[MPa]`"
    b1::FT = 1
    "C of first part"
    c1::FT = 5
    "F of first part"
    f1::FT = 0.5
    "B of second part `[MPa]`"
    b2::FT = 2
    "C of second part"
    c2::FT = 5
    "F of second part"
    f2::FT = 1 - f1
end




"""
    $(TYPEDEF)

Struct that contains single Weibull function parameters.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct WeibullSingle{FT} <: AbstractXylemVC{FT}
    "B of first part `[MPa]`"
    b::FT = 2
    "C of first part"
    c::FT = 5
end
=#
