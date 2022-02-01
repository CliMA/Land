###############################################################################
#
# Flow mode, steady state or non-steady state
#
###############################################################################
"""
    abstract type AbstractFlowMode

Hierarchy of AbstractFlowMode
- [`SteadyStateMode`](@ref)
- [`NonSteadyStateMode`](@ref)
"""
abstract type AbstractFlowMode end




"""
    struct SteadyStateMode <: AbstractFlowMode
"""
struct SteadyStateMode <: AbstractFlowMode end




"""
    struct NonSteadyStateMode <: AbstractFlowMode
"""
struct NonSteadyStateMode <: AbstractFlowMode end
