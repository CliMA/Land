###############################################################################
#
# Photosynthesis mode
#
###############################################################################
"""
    abstract type AbstractCalculationMode

Hierarchy of AbstractCalculationMode
- [`GCO₂Mode`](@ref)
- [`PCO₂Mode`](@ref)
"""
abstract type AbstractCalculationMode end




"""
    struct GCO₂Mode <: AbstractCalculationMode

Calculate leaf photosynthesis using leaf internal CO₂ partial pressure
"""
struct GCO₂Mode <: AbstractCalculationMode end




"""
    struct PCO₂Mode <: AbstractCalculationMode

Calculate leaf photosynthesis using total leaf diffusive conductance to CO₂
"""
struct PCO₂Mode <: AbstractCalculationMode end
