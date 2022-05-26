#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-May-25: add abstract type for soil-plant-air continuum
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractSPACSystem:
"""
abstract type AbstractSPACSystem{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-25: toy SPAC system
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for simplest SPAC system

# Fields

$(TYPEDFIELDS)

"""
mutable struct MonoElementSAPC{FT} <: AbstractSPACSystem{FT}
    # parameters that do not change with time
    "Leaf hydrualic system"
    LEAF::Leaf{FT}
    "Root hydraulic system"
    ROOT::RootHydraulics{FT}
    "Stem hydraulic system"
    STEM::StemHydraulics{FT}

    # caches to speed up calculations
    "Relative hydraulic conductance"
    _krs::Array{FT,1}
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-25: SPAC system for monospecies grass
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for monospecies grass SPAC system

# Fields

$(TYPEDFIELDS)

"""
mutable struct MonoGrassSAPC{FT} <: AbstractSPACSystem{FT}
    # parameters that do not change with time
    "Leaf hydrualic system"
    LEAVES::Vector{Leaf{FT}}
    "Corresponding air layer per canopy layer"
    LEAVES_INDEX::Vector{Int}
    "Number of canopy layers"
    N_CANOPY::Int
    "Number of root layers"
    N_ROOT::Int
    "Root hydraulic system"
    ROOTS::Vector{RootHydraulics{FT}}
    "Corresponding soil layer per root layer"
    ROOTS_INDEX::Vector{Int}

    # caches to speed up calculations
    "Conductances for each root layer at given flow"
    _ks::Array{FT,1}
    "Pressure for each root layer at given flow"
    _ps::Array{FT,1}
    "Flow rate per root layer"
    _qs::Array{FT,1}
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-25: SPAC system for monospecies palm
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for monospecies palm SPAC system (with trunk)

# Fields

$(TYPEDFIELDS)

"""
mutable struct MonoPalmSAPC{FT} <: AbstractSPACSystem{FT}
    # parameters that do not change with time
    "Leaf hydrualic system"
    LEAVES::Vector{Leaf{FT}}
    "Corresponding air layer per canopy layer"
    LEAVES_INDEX::Vector{Int}
    "Number of canopy layers"
    N_CANOPY::Int
    "Number of root layers"
    N_ROOT::Int
    "Root hydraulic system"
    ROOTS::Vector{RootHydraulics{FT}}
    "Corresponding soil layer per root layer"
    ROOTS_INDEX::Vector{Int}
    "Trunk hydraulic system"
    TRUNK::StemHydraulics{FT}

    # caches to speed up calculations
    "Conductances for each root layer at given flow"
    _ks::Array{FT,1}
    "Pressure for each root layer at given flow"
    _ps::Array{FT,1}
    "Flow rate per root layer"
    _qs::Array{FT,1}
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-25: SPAC system for monospecies tree
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for monospecies tree SPAC system (with trunk and branches)

# Fields

$(TYPEDFIELDS)

"""
mutable struct MonoTreeSAPC{FT} <: AbstractSPACSystem{FT}
    # parameters that do not change with time
    "Branch hydraulic system"
    BRANCHES::Vector{StemHydraulics{FT}}
    "Leaf hydrualic system"
    LEAVES::Vector{Leaf{FT}}
    "Corresponding air layer per canopy layer"
    LEAVES_INDEX::Vector{Int}
    "Number of canopy layers"
    N_CANOPY::Int
    "Number of root layers"
    N_ROOT::Int
    "Root hydraulic system"
    ROOTS::Vector{RootHydraulics{FT}}
    "Corresponding soil layer per root layer"
    ROOTS_INDEX::Vector{Int}
    "Trunk hydraulic system"
    TRUNK::StemHydraulics{FT}

    # caches to speed up calculations
    "Conductances for each root layer at given flow"
    _ks::Array{FT,1}
    "Pressure for each root layer at given flow"
    _ps::Array{FT,1}
    "Flow rate per root layer"
    _qs::Array{FT,1}
end
