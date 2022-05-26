#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-May-25: add Root structure
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save stem parameters

# Fields

$(TYPEDFIELDS)

"""
mutable struct Stem{FT<:AbstractFloat}
    # parameters that do not change with time
    "[`StemHydraulics`](@ref) type stem hydraulic system"
    HS::StemHydraulics{FT}

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

    Stem{FT}() where {FT<:AbstractFloat}

Construct a Stem structure
"""
Stem{FT}() where {FT<:AbstractFloat} = Stem{FT}(StemHydraulics{FT}(), T_25());
