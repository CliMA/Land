###############################################################################
#
# Leaf boundary layer parameter set
#
###############################################################################
#= AbstractLeafBLParaSet type tree
---> LeafBLParaSetFixed
---> LeafBLParaSetGentine
=#
"""
    AbstractLeafBLParaSet

Hierarchy of the `AbstractLeafBLParaSet`:
- [`LeafBLParaSetFixed`](@ref)
- [`LeafBLParaSetGentine`](@ref)
"""
abstract type AbstractLeafBLParaSet end




"""
    struct LeafBLParaSetFixed{FT<:AbstractFloat}

Leaf boundary layer with fixed resistance

# Fields
$(DocStringExtensions.FIELDS)
"""
struct LeafBLParaSetFixed{FT<:AbstractFloat} <: AbstractLeafBLParaSet
    "Resistance of the leaf boundary layer `[m² s mol⁻¹]`"
    ra::FT
end




"""
    struct LeafBLParaSetGentine

Gentine's boundary layer scheme that computes Monin Obhukow length and
    resistances across the leaf to Ca scale
"""
struct LeafBLParaSetGentine <: AbstractLeafBLParaSet end
