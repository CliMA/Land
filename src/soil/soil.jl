#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-08: add Soil structure
#     2022-Jun-09: add fields: e_net_diffuse, e_net_direct
#     2022-Jun-10: add fields: r_net_sw
# To do:
#     TODO: abstractize the albedo types (including leaf)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for Soil

# Fields

$(TYPEDFIELDS)

"""
mutable struct Soil{FT<:AbstractFloat}
    # parameters that do not change with time
    "Soil moisture retention curve"
    VC::AbstractSoilVC{FT}
    "Mean depth"
    Z::FT
    "Depth boundaries"
    ZS::Vector{FT}

    # prognostic variables that change with time
    "Temperature"
    t::FT
    "Soil water content"
    θ::FT

    # diagnostic variables that change with time
    "Net diffuse radiation at each canopy layer `[mW m⁻² nm⁻¹]`"
    e_net_diffuse::Vector{FT}
    "Net direct radiation at each canopy layer `[mW m⁻² nm⁻¹]`"
    e_net_direct::Vector{FT}
    "Net shortwave energy absorption `[W m⁻²]`"
    r_net_sw::FT
    "Reflectance"
    ρ_sw::Vector{FT}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-08: add constructor
#     2022-Jun-09: add fields: e_net_diffuse, e_net_direct
#     2022-Jun-10: add fields: r_net_sw
#
#######################################################################################################################################################################################################
"""

    Soil{FT}(zs::Vector{FT}; soil_type::String = "Loam", n_λ::Int = 114) where {FT<:AbstractFloat}

Construct a soil struct, given
- `zs` Soil upper and lower boundaries
- `soil_type` Soil type name
- `n_λ` Number of shortwave wavelength bins
"""
Soil{FT}(zs::Vector{FT}; soil_type::String = "Loam", n_λ::Int = 114) where {FT<:AbstractFloat} = (
    _svc = VanGenuchten{FT}(soil_type);
    return Soil{FT}(
                _svc,           # VC
                mean(zs),       # Z
                zs,             # ZS
                T_25(FT),       # t
                _svc.Θ_SAT,     # θ
                zeros(FT,n_λ),  # e_net_diffuse
                zeros(FT,n_λ),  # e_net_direct
                0,              # r_net_sw
                zeros(FT,n_λ)   # ρ_sw
    )
);
