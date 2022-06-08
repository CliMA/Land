#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-08: add Soil structure
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
    "Reflectance"
    ρ_sw::Vector{FT}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-08: add constructor
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
                zeros(FT,n_λ)   # ρ_sw
    )
);
