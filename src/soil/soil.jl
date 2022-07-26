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
#     2022-Jun-14: make soil albedo a two-element vector for PAR and NIR
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for broadband soil albedo

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct BroadbandSoilAlbedo{FT<:AbstractFloat} <: AbstractSoilAlbedo{FT}
    # Constants
    "Reflectance for longwave radiation"
    ρ_LW::FT = 0.06

    # Diagnostic variables
    "Net diffuse radiation at top soil `[W m⁻²]`"
    e_net_diffuse::FT = 0
    "Net direct radiation at top soil `[W m⁻²]`"
    e_net_direct::FT = 0
    "Net longwave energy absorption `[W m⁻²]`"
    r_net_lw::FT = 0
    "Net shortwave energy absorption `[W m⁻²]`"
    r_net_sw::FT = 0
    "Reflectance for shortwave radiation (for PAR and NIR)"
    ρ_sw::Vector{FT} = FT[0, 0]
end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-14: add struct for hyperspectral soil albedo
#     2022-Jun-14: add constructor
#     2022-Jun-14: add fields to compute soil hyperspectral albedo in CanopyRadiativeTransfer.jl
#     2022-Jun-14: add wls in constructor function and remove n_λ
#     2022-Jul-22: rename Ρ (greek) to ρ
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for hyperspectral soil albedo

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct HyperspectralSoilAlbedo{FT<:AbstractFloat} <: AbstractSoilAlbedo{FT}
    # File path to the Netcdf dataset
    "File path to the Netcdf dataset"
    DATASET::String = LAND_2021

    # Dimensions
    "Number of wavelength bins for NIR"
    DIM_NIR::Int = 79
    "Number of wavelength bins"
    DIM_WL::Int = 114

    # Constants
    "A matrix of characteristic curves"
    MAT_ρ::Matrix{FT} = FT[read_nc(DATASET, "GSV_1") read_nc(DATASET, "GSV_2") read_nc(DATASET, "GSV_3") read_nc(DATASET, "GSV_4")]
    "Reflectance for longwave radiation"
    ρ_LW::FT = 0.06

    # Diagnostic variables
    "Net diffuse radiation at top soil `[mW m⁻² nm⁻¹]`"
    e_net_diffuse::Vector{FT} = zeros(FT, DIM_WL)
    "Net direct radiation at top soil `[mW m⁻² nm⁻¹]`"
    e_net_direct::Vector{FT} = zeros(FT, DIM_WL)
    "Net longwave energy absorption `[W m⁻²]`"
    r_net_lw::FT = 0
    "Net shortwave energy absorption `[W m⁻²]`"
    r_net_sw::FT = 0
    "Reflectance for shortwave radiation"
    ρ_sw::Vector{FT} = zeros(FT, DIM_WL)

    # Cache variables
    "Cache variable with length of NIR"
    _tmp_vec_nir::Vector{FT} = zeros(FT, DIM_NIR)
    "Weights of the four characteristic curves"
    _weight::Vector{FT} = zeros(FT, 4)
    "Cache variable to store ρ_PAR and ρ_NIR (a segmented curve)"
    _ρ_sw::Vector{FT} = zeros(FT, DIM_WL)
end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jul-13: add SoilLayer structure
#     2022-Jul-13: add field K_MAX, K_REF, k, ψ, and ∂θ∂t
#     2022-Jul-14: remove field K_REF
#     2022-Jul-14: add field ∂G∂t (renamed to ∂e∂t), ΔZ
#     2022-Jul-26: move field K_MAX to structure VanGenuchten and BrooksCorey
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for soil layer

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SoilLayer{FT<:AbstractFloat}
    # Constants
    "Specific heat capacity of soil `[J K⁻¹ kg⁻¹]`"
    CP::FT = 760
    "Soil thermal conductivity `[W m⁻¹ K⁻¹]`"
    Λ_THERMAL::FT = 3
    "Dry soil density `[kg m⁻³]`"
    ρ::FT = 2650

    # Embedded structures
    "Soil moisture retention curve"
    VC::Union{BrooksCorey{FT}, VanGenuchten{FT}} = VanGenuchten{FT}("Loam")

    # Geometry information
    "Depth boundaries `[m]`"
    ZS::Vector{FT} = FT[0,-1]
    "Mean depth `[m]`"
    Z::FT = (ZS[1] + ZS[2]) / 2
    "Layer thickness `[m]`"
    ΔZ::FT = ZS[1] - ZS[2]

    # Prognostic variables (not used for ∂y∂t)
    "Temperature `[K]`"
    t::FT = T₂₅()

    # Prognostic variables (used for ∂y∂t)
    "Total stored energy per volume `[J m⁻³]`"
    e::FT = (CP * ρ + VC.Θ_SAT * CP_L() * ρ_H₂O()) * t
    "Soil water content"
    θ::FT = VC.Θ_SAT
    "Marginal increase in energy `[W m⁻²]`"
    ∂e∂t::FT = 0
    "Marginal increase in soil water content `[s⁻¹]`"
    ∂θ∂t::FT = 0

    # Diagnostic variables
    "Soil hydraulic conductance per area `[mol m⁻² s⁻¹ MPa⁻¹]`"
    k::FT = 0
    "Matric potential `[MPa]`"
    ψ::FT = 0

    # Cache variables
    "Combined specific heat capacity of soil `[J K⁻¹ kg⁻¹]`"
    _cp::FT = 0
    "Combined soil thermal conductance `[W m⁻² K⁻¹]`"
    _λ_thermal::FT = 0
end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-08: add Soil structure
#     2022-Jun-08: add constructor
#     2022-Jun-09: add fields: e_net_diffuse, e_net_direct
#     2022-Jun-10: add fields: r_net_lw, r_net_sw, ρ_lw
#     2022-Jun-14: add abstractized soil albedo
#     2022-Jun-13: use Union instead of Abstract... for type definition
#     2022-Jun-14: add field for soil color class
#     2022-Jul-13: move VC, Z, t, and θ to SoilLayer
#     2022-Jul-13: add field AREA, _k, _q, and _δψ
#     2022-Jul-13: add field _λ_thermal, _q_thermal, and _δt
#     2022-Jul-15: add field runoff
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for Soil

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Soil{FT<:AbstractFloat}
    # Dimensions
    "Dimension of soil layers"
    DIM_SOIL::Int = 5

    # General information
    "Total area of the soil `[m²]`"
    AREA::FT = 100
    "Color class as in CLM"
    COLOR::Int = 1
    "Soil layers boundaries"
    ZS::Vector{FT} = FT[0, -0.1, -0.2, -0.5, -1, -2]

    # Embedded structures
    "Albedo related structure"
    ALBEDO::Union{BroadbandSoilAlbedo{FT}, HyperspectralSoilAlbedo{FT}} = HyperspectralSoilAlbedo{FT}()
    "Soil layers"
    LAYERS::Vector{SoilLayer{FT}} = SoilLayer{FT}[SoilLayer{FT}(VC = VanGenuchten{FT}("Loam"), ZS = ZS[_i:_i+1]) for _i in 1:length(ZS)-1]

    # Diagnostic variables
    "Surface runoff due to heavy precipitation during the time step `[mol m⁻²]`"
    runoff::FT = 0

    # Cache variables
    "Soil hydraulic conductance between layers per area `[mol m⁻² s⁻¹ MPa⁻¹]`"
    _k::Vector{FT} = zeros(FT, DIM_SOIL - 1)
    "Flux between layers per area `[mol m⁻² s⁻¹]`"
    _q::Vector{FT} = zeros(FT, DIM_SOIL - 1)
    "Thermal flux between layers per area `[mol m⁻² s⁻¹]`"
    _q_thermal::Vector{FT} = zeros(FT, DIM_SOIL - 1)
    "Soil temperature difference between layers `[MPa]`"
    _δt::Vector{FT} = zeros(FT, DIM_SOIL - 1)
    "Soil metric potential difference between layers `[MPa]`"
    _δψ::Vector{FT} = zeros(FT, DIM_SOIL - 1)
    "Soil thermal conductance between layers per area `[W m⁻² K⁻¹]`"
    _λ_thermal::Vector{FT} = zeros(FT, DIM_SOIL - 1)
end
