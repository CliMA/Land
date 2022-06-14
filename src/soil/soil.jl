#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jun-14: add abstract type for soil albedo
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractSoilAlbedo:
- [`BroadbandSoilAlbedo`](@ref)
- [`HyperspectralSoilAlbedo`](@ref)
"""
abstract type AbstractSoilAlbedo{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-14: add struct for broadband soil albedo
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for broadband soil albedo

# Fields

$(TYPEDFIELDS)

"""
mutable struct BroadbandSoilAlbedo{FT} <: AbstractSoilAlbedo{FT}
    # diagnostic variables that change with time
    "Net diffuse radiation at top soil `[W m⁻²]`"
    e_net_diffuse::FT
    "Net direct radiation at top soil `[W m⁻²]`"
    e_net_direct::FT
    "Net longwave energy absorption `[W m⁻²]`"
    r_net_lw::FT
    "Net shortwave energy absorption `[W m⁻²]`"
    r_net_sw::FT
    "Reflectance for longwave radiation"
    ρ_lw::FT
    "Reflectance for shortwave radiation"
    ρ_sw::FT
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-14: add constructor
#
#######################################################################################################################################################################################################
"""

    BroadbandSoilAlbedo{FT}() where {FT<:AbstractFloat}

Construct a broadband soil albedo struct
"""
BroadbandSoilAlbedo{FT}() where {FT<:AbstractFloat} = (
    return BroadbandSoilAlbedo{FT}(
                0,      # e_net_diffuse
                0,      # e_net_direct
                0,      # r_net_lw
                0,      # r_net_sw
                0.06,   # ρ_lw
                0       # ρ_sw
    )
);


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-14: add struct for hyperspectral soil albedo
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for hyperspectral soil albedo

# Fields

$(TYPEDFIELDS)

"""
mutable struct HyperspectralSoilAlbedo{FT} <: AbstractSoilAlbedo{FT}
    # diagnostic variables that change with time
    "Net diffuse radiation at top soil `[mW m⁻² nm⁻¹]`"
    e_net_diffuse::Vector{FT}
    "Net direct radiation at top soil `[mW m⁻² nm⁻¹]`"
    e_net_direct::Vector{FT}
    "Net longwave energy absorption `[W m⁻²]`"
    r_net_lw::FT
    "Net shortwave energy absorption `[W m⁻²]`"
    r_net_sw::FT
    "Reflectance for longwave radiation"
    ρ_lw::FT
    "Reflectance for shortwave radiation"
    ρ_sw::Vector{FT}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-14: add constructor
#
#######################################################################################################################################################################################################
"""

    HyperspectralSoilAlbedo{FT}(; n_λ::Int = 114) where {FT<:AbstractFloat}

Construct a hyperspectral soil albedo struct, given
- `n_λ` Number of shortwave wavelength bins
"""
HyperspectralSoilAlbedo{FT}(; n_λ::Int = 114) where {FT<:AbstractFloat} = (
    return HyperspectralSoilAlbedo{FT}(
                zeros(FT,n_λ),  # e_net_diffuse
                zeros(FT,n_λ),  # e_net_direct
                0,              # r_net_lw
                0,              # r_net_sw
                0.06,           # ρ_lw
                zeros(FT,n_λ)   # ρ_sw
    )
);


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-08: add Soil structure
#     2022-Jun-09: add fields: e_net_diffuse, e_net_direct
#     2022-Jun-10: add fields: r_net_lw, r_net_sw, ρ_lw
#     2022-Jun-14: add abstractized soil albedo
#     2022-Jun-13: use Union instead of Abstract... for type definition
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
    "Albedo related structure"
    ALBEDO::Union{BroadbandSoilAlbedo{FT}, HyperspectralSoilAlbedo{FT}}
    "Soil moisture retention curve"
    VC::Union{BrooksCorey{FT}, VanGenuchten{FT}}
    "Mean depth"
    Z::FT
    "Depth boundaries"
    ZS::Vector{FT}

    # prognostic variables that change with time
    "Temperature"
    t::FT
    "Soil water content"
    θ::FT
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-08: add constructor
#     2022-Jun-09: add fields: e_net_diffuse, e_net_direct
#     2022-Jun-10: add fields: r_net_lw, r_net_sw, ρ_lw
#     2022-Jun-14: add abstractized soil albedo
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
    _sab = n_λ > 1 ? HyperspectralSoilAlbedo{FT}(n_λ = n_λ) : BroadbandSoilAlbedo{FT}();

    return Soil{FT}(
                _sab,       # ALBEDO
                _svc,       # VC
                mean(zs),   # Z
                zs,         # ZS
                T_25(FT),   # t
                _svc.Θ_SAT  # θ
    )
);
