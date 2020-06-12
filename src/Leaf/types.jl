###############################################################################
#
# Leaf boundary layer parameter set
# These types are documented in the Leaf page
#
###############################################################################
"""
    AbstractLeafBLParaSet

A parameter set that stores boundary layer parameters
"""
abstract type AbstractLeafBLParaSet end




"""
    struct LeafBLParaSetFixed{FT, resistance}

A leaf boundary layer with fixed resistance, namely ensuring surface concentration
"""
Base.@kwdef struct LeafBLParaSetFixed{FT, resistance} <: AbstractLeafBLParaSet
    ra::FT = resistance
end




"""
    struct LeafBLParaSetGentine

Gentine's boundary layer scheme that computes Monin Obhukow length and resistances across the leaf to Ca scale
"""
struct LeafBLParaSetGentine <: AbstractLeafBLParaSet end
