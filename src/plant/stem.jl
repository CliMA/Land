#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-May-25: add Stem structure
#     2022-Jul-15: add fields e, ∂e∂t
#     2022-Jul-19: use kwdef for the constructor
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save stem parameters

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Stem{FT<:AbstractFloat}
    # parameters that do not change with time
    "[`StemHydraulics`](@ref) type stem hydraulic system"
    HS::StemHydraulics{FT} = StemHydraulics{FT}()

    # prognostic variables that change with time (# TODO: add wood storage as well)
    "Total stored energy in water `[J]`"
    e::FT = T_25() * sum(HS.v_storage) * CP_L_MOL(FT)
    "Current temperature"
    t::FT = T_25()

    # diagnostic variables that change with time
    "Marginal increase in energy `[W]`"
    ∂e∂t::FT = 0
end
