###############################################################################
#
# Plant Hydraulic system
#
###############################################################################
"""
    abstract type AbstractPlantOrganism{FT}

Hierachy of AbstractPlantOrganism
- [`GrassLikeOrganism`](@ref)
- [`PalmLikeOrganism`](@ref)
- [`TreeLikeOrganism`](@ref)
"""
abstract type AbstractPlantOrganism{FT<:AbstractFloat} end




"""
    mutable struct GrassLikeOrganism{FT}

A plant hydraulic system like a grass, which contains multiple root layers, and
multiple canopy layers. No trunk or branch system applies.

# Fields
$(TYPEDFIELDS)
"""
mutable struct GrassLikeOrganism{FT} <: AbstractPlantOrganism{FT}
    # structre information
    "Root Layers"
    n_root  ::Int
    "Canopy Layers"
    n_canopy::Int

    # Arrays of roots and leaves
    "Roots system"
    roots ::Array{RootHydraulics{FT},1}
    "Leaves"
    leaves::Array{LeafHydraulics{FT},1}

    # Root and canopy index in Soil and Atmosphere
    "Corresponding soil layer per root layer"
    root_index_in_soil ::Array{Int,1}
    "Corresponding air layer per canopy layer"
    canopy_index_in_air::Array{Int,1}

    # containers for root flows
    "Conductances for each layer at given flow"
    cache_k::Array{FT,1}
    "Pressure for each layer at given flow"
    cache_p::Array{FT,1}
    "Flow rate"
    cache_q::Array{FT,1}
end




"""
    mutable struct PalmLikeOrganism{FT}

A plant hydraulic system like a palm, which contains multiple root layers, one
trunk, and multiple canopy layers. No branch system applies.

# Fields
$(TYPEDFIELDS)
"""
mutable struct PalmLikeOrganism{FT} <: AbstractPlantOrganism{FT}
    # structre information
    "Root Layers"
    n_root  ::Int
    "Canopy Layers"
    n_canopy::Int

    # Arrays of roots and leaves
    "Roots system"
    roots ::Array{RootHydraulics{FT},1}
    "Trunk"
    trunk ::StemHydraulics{FT}
    "Leaves"
    leaves::Array{LeafHydraulics{FT},1}

    # Root and canopy index in Soil and Atmosphere
    "Corresponding soil layer per root layer"
    root_index_in_soil ::Array{Int,1}
    "Corresponding air layer per canopy layer"
    canopy_index_in_air::Array{Int,1}

    # containers for root flows
    "Conductances for each layer at given flow"
    cache_k::Array{FT,1}
    "Pressure for each layer at given flow"
    cache_p::Array{FT,1}
    "Flow rate"
    cache_q::Array{FT,1}
end




"""
    mutable struct TreeLikeOrganism{FT}

A plant hydraulic system like a tree, which contains multiple root layers, one
trunk, and multiple branch and canopy layers.

# Fields
$(TYPEDFIELDS)
"""
mutable struct TreeLikeOrganism{FT} <: AbstractPlantOrganism{FT}
    # structre information
    "Root Layers"
    n_root  ::Int
    "Canopy Layers"
    n_canopy::Int

    # Arrays of roots and leaves
    "Roots system"
    roots ::Array{RootHydraulics{FT},1}
    "Trunk"
    trunk ::StemHydraulics{FT}
    "Branch system"
    branch::Array{StemHydraulics{FT},1}
    "Leaves"
    leaves::Array{LeafHydraulics{FT},1}

    # Root and canopy index in Soil and Atmosphere
    "Corresponding soil layer per root layer"
    root_index_in_soil ::Array{Int,1}
    "Corresponding air layer per canopy layer"
    canopy_index_in_air::Array{Int,1}

    # containers for root flows
    "Conductances for each layer at given flow"
    cache_k::Array{FT,1}
    "Pressure for each layer at given flow"
    cache_p::Array{FT,1}
    "Flow rate"
    cache_q::Array{FT,1}
end




"""
    mutable struct TreeSimple{FT}

A plant hydraulic system with one root, one stem, and one leaf for testing
    purpose

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct TreeSimple{FT} <: AbstractPlantOrganism{FT}
    # Arrays of roots and leaves
    "Root"
    root::RootHydraulics{FT} = RootHydraulics{FT}()
    "Stem"
    stem::StemHydraulics{FT} = StemHydraulics{FT}()
    "Leaf"
    leaf::LeafHydraulics{FT} = LeafHydraulics{FT}()

    # Local container for tree information
    "Relative hydraulic conductance"
    krs::Array{FT,1} = ones(FT,4)
end
