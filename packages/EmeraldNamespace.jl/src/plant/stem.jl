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
    # Embedded structures
    "[`StemHydraulics`](@ref) type stem hydraulic system"
    HS::StemHydraulics{FT} = StemHydraulics{FT}()

    # Prognostic variables (not used for ∂y∂t)
    "Current temperature"
    t::FT = T₂₅()

    # Prognostic variables (used for ∂y∂t)
    "Total stored energy in water `[J]`" # TODO: add wood storage as well
    e::FT = T₂₅() * sum(HS.v_storage) * CP_L_MOL(FT)
    "Marginal increase in energy `[W]`"
    ∂e∂t::FT = 0
end
