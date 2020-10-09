###############################################################################
#
# Leaf fluorescence-related parameter set
#
###############################################################################
#= AbstractFluoModelParaSet type tree
---> FluoParaSet
=#
"""
    abstract type AbstractFluoModelParaSet{FT}

Hierarchy of the `AbstractFluoModelParaSet`:
- [`FluoParaSet`](@ref)
"""
abstract type AbstractFluoModelParaSet{FT} end




"""
    mutable struct FluoParaSet{FT}

A `AbstractFluoModelParaSet` type paramter set.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct FluoParaSet{FT<:AbstractFloat} <: AbstractFluoModelParaSet{FT}
    "Fluorescence model coefficient"
    Kn1::FT
    "Fluorescence model coefficient"
    Kn2::FT
    "Fluorescence model coefficient"
    Kn3::FT
end
