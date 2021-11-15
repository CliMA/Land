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
    mutable struct CytoFluoParaSet{FT}

A `AbstractFluoModelParaSet` type paramter set.

# Fields
$(TYPEDFIELDS)
"""
struct CytoFluoParaSet{FT<:AbstractFloat} <: AbstractFluoModelParaSet{FT}
    "Fluorescence model coefficient"
    Kr1::FT
    "Fluorescence model coefficient"
    Kr2::FT
    "Fluorescence model coefficient"
    Kr3::FT
end




"""
    mutable struct FluoParaSet{FT}

A `AbstractFluoModelParaSet` type paramter set.

# Fields
$(TYPEDFIELDS)
"""
struct FluoParaSet{FT<:AbstractFloat} <: AbstractFluoModelParaSet{FT}
    "Fluorescence model coefficient"
    Kr1::FT
    "Fluorescence model coefficient"
    Kr2::FT
    "Fluorescence model coefficient"
    Kr3::FT
end
