#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-May-25: add Root structure
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save root parameters

# Fields

$(TYPEDFIELDS)

"""
mutable struct Root{FT<:AbstractFloat}
    # parameters that do not change with time
    "[`RootHydraulics`](@ref) type root hydraulic system"
    HS::RootHydraulics{FT}

    # prognostic variables that change with time
    "Current temperature"
    t::FT
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-May-25: add constructor
#
#######################################################################################################################################################################################################
"""

    Root{FT}() where {FT<:AbstractFloat}

Construct a Root structure
"""
Root{FT}() where {FT<:AbstractFloat} = Root{FT}(RootHydraulics{FT}(), T_25());
