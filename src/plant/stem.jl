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

    Stem{FT}(; ssm::Bool = true) where {FT<:AbstractFloat}

Construct a Stem structure, given
- `ssm` Whether the flow rate is at steady state
"""
Stem{FT}(; ssm::Bool = true) where {FT<:AbstractFloat} = (
    _hs = StemHydraulics{FT}(ssm = ssm);

    return Stem{FT}(
                _hs,                    # HS
                T_25() * _hs.v_storage, # e
                T_25(),                 # t
                0                       # ∂e∂t
    )
);
