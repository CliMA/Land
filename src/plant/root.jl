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

    # prognostic variables that change with time (# TODO: add wood storage as well)
    "Total stored energy in water `[J]`"
    e::FT
    "Current temperature"
    t::FT

    # diagnostic variables that change with time
    "Marginal increase in energy `[W]`"
    ∂e∂t::FT
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-May-25: add constructor
#     2022-May-31: add steady state mode option to input options
#
#######################################################################################################################################################################################################
"""

    Root{FT}(; ssm::Bool = true) where {FT<:AbstractFloat}

Construct a Root structure, given
- `ssm` Whether the flow rate is at steady state
"""
Root{FT}(; ssm::Bool = true) where {FT<:AbstractFloat} = (
    _hs = RootHydraulics{FT}(ssm = ssm);

    return Root{FT}(
                _hs,                    # HS
                T_25() * _hs.v_storage, # e
                T_25(),                 # t
                0                       # ∂e∂t
    )
);
