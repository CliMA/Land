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
#     2022-Jul-13: use @kwdef for the constructor
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for broadband soil albedo

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct BroadbandSoilAlbedo{FT<:AbstractFloat} <: AbstractSoilAlbedo{FT}
    # diagnostic variables that change with time
    "Net diffuse radiation at top soil `[W m⁻²]`"
    e_net_diffuse::FT = 0
    "Net direct radiation at top soil `[W m⁻²]`"
    e_net_direct::FT = 0
    "Net longwave energy absorption `[W m⁻²]`"
    r_net_lw::FT = 0
    "Net shortwave energy absorption `[W m⁻²]`"
    r_net_sw::FT = 0
    "Reflectance for longwave radiation"
    ρ_lw::FT = 0.06
    "Reflectance for shortwave radiation (for PAR and NIR)"
    ρ_sw::Vector{FT} = FT[0,0]
end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-14: add struct for hyperspectral soil albedo
#     2022-Jun-14: add constructor
#     2022-Jun-14: add fields to compute soil hyperspectral albedo in CanopyRadiativeTransfer.jl
#     2022-Jun-14: add wls in constructor function and remove n_λ
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

    # dimensions
    "Number of wavelength bins for NIR"
    DIM_NIR::Int = 79
    "Number of wavelength bins"
    DIM_WL::Int = 114

    # parameters that do not change with time
    "A matrix of characteristic curves"
    MAT_Ρ::Matrix{FT} = FT[read_nc(DATASET, "GSV_1") read_nc(DATASET, "GSV_2") read_nc(DATASET, "GSV_3") read_nc(DATASET, "GSV_4")]
    "Reflectance for longwave radiation"
    Ρ_LW::FT = 0.06

    # diagnostic variables that change with time
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

    # caches to speed up calculations
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
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for soil layer

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SoilLayer{FT<:AbstractFloat}
    # parameters that do not change with time
    "Specific heat capacity of soil `[J K⁻¹ kg⁻¹]`"
    CP::FT = 760
    "Maximum soil hydraulic conductivity at 25 °C `[mol m⁻¹ s⁻¹ MPa⁻¹]`"
    K_MAX::FT = 10000
    "Soil moisture retention curve"
    VC::Union{BrooksCorey{FT}, VanGenuchten{FT}} = VanGenuchten{FT}("Loam")
    "Mean depth"
    Z::FT = -0.5
    "Depth boundaries"
    ZS::Vector{FT} = FT[0,-1]
    "Layer thickness `[m]`"
    ΔZ::FT = ZS[1] - ZS[2]
    "Dry soil density `[kg m⁻³]`"
    Ρ::FT = 2650
    "Soil thermal conductivity `[W m⁻¹ K⁻¹]`"
    Λ_THERMAL::FT = 3

    # prognostic variables that change with time
    "Total stored energy per volume `[J m⁻³]`"
    e::FT = T_25() * (CP * Ρ + VC.Θ_SAT * CP_L() * ρ_H₂O())
    "Soil water content"
    θ::FT = VC.Θ_SAT

    # diagnostic variables that change with time
    "Combined specific heat capacity of soil `[J K⁻¹ kg⁻¹]`"
    cp::FT = 0
    "Soil hydraulic conductance per area `[mol m⁻² s⁻¹ MPa⁻¹]`"
    k::FT = 0
    "Temperature `[K]`"
    t::FT = T_25()
    "Combined soil thermal conductance `[W m⁻² K⁻¹]`"
    λ_thermal::FT = 0
    "Matric potential `[MPa]`"
    ψ::FT = 0
    "Marginal increase in energy `[W m⁻²]`"
    ∂e∂t::FT = 0
    "Marginal increase in soil water content `[s⁻¹]`"
    ∂θ∂t::FT = 0
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
#     2022-Jun-14: separate the constructor for hyperspectral albedo
#     2022-Jun-14: separate the constructor for broadband albedo
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
mutable struct Soil{FT<:AbstractFloat}
    # parameters that do not change with time
    "Albedo related structure"
    ALBEDO::Union{BroadbandSoilAlbedo{FT}, HyperspectralSoilAlbedo{FT}}
    "Total area of the soil `[m²]`"
    AREA::FT
    "Color class as in CLM"
    COLOR::Int
    "Soil layers"
    LAYERS::Vector{SoilLayer{FT}}
    "Number of soil layers"
    N_LAYER::Int

    # diagnostic parameters that changes with time
    "Surface runoff due to heavy precipitation during the time step `[mol m⁻²]`"
    runoff::FT

    # cache used for calculations
    "Soil hydraulic conductance between layers per area `[mol m⁻² s⁻¹ MPa⁻¹]`"
    _k::Vector{FT}
    "Flux between layers per area `[mol m⁻² s⁻¹]`"
    _q::Vector{FT}
    "Thermal flux between layers per area `[mol m⁻² s⁻¹]`"
    _q_thermal::Vector{FT}
    "Soil temperature difference between layers `[MPa]`"
    _δt::Vector{FT}
    "Soil metric potential difference between layers `[MPa]`"
    _δψ::Vector{FT}
    "Soil thermal conductance between layers per area `[W m⁻² K⁻¹]`"
    _λ_thermal::Vector{FT}
end


"""

    Soil{FT}(zs::Vector, area::Number, wls::WaveLengthSet{FT} = WaveLengthSet{FT}(); soil_type::String = "Loam") where {FT<:AbstractFloat}

Construct a soil struct with hyperspectral albedo, given
- `zs` Soil upper and lower boundaries
- `area` Surface area of the soil (per tree for MonoML*SPAC)
- `wls` [`WaveLengthSet`](@ref) type struct that defines wavelength settings
- `soil_type` Soil type name

"""
Soil{FT}(zs::Vector, area::Number, wls::WaveLengthSet{FT} = WaveLengthSet{FT}(); soil_type::String = "Loam") where {FT<:AbstractFloat} = (
    _soil = Soil{FT}(zs, area, true; soil_type = soil_type);
    _soil.ALBEDO = HyperspectralSoilAlbedo{FT}();

    return _soil
);


"""

    Soil{FT}(zs::Vector, area::Number, broadband::Bool; soil_type::String = "Loam") where {FT<:AbstractFloat}

Construct a soil struct with broadband albedo, given
- `zs` Soil upper and lower boundaries
- `area` Surface area of the soil (per tree for MonoML*SPAC)
- `broadband` Indicating broadband soil albedo
- `soil_type` Soil type name

"""
Soil{FT}(zs::Vector, area::Number, broadband::Bool; soil_type::String = "Loam") where {FT<:AbstractFloat} = (
    _layers = SoilLayer{FT}[];
    _n_layer = length(zs) - 1;
    for _i in 1:_n_layer
        push!(_layers, SoilLayer{FT}(K_MAX = 1e4, VC = VanGenuchten{FT}(soil_type), Z = mean(zs[_i:_i+1]), ZS = zs[_i:_i+1]));
    end;

    _sab = BroadbandSoilAlbedo{FT}();

    return Soil{FT}(
                _sab,               # ALBEDO
                100,                # AREA
                1,                  # COLOR
                _layers,            # LAYERS
                _n_layer,           # N_LAYER
                0,                  # runoff
                zeros(_n_layer-1),  # _k
                zeros(_n_layer-1),  # _q
                zeros(_n_layer-1),  # _q_thermal
                zeros(_n_layer-1),  # _δt
                zeros(_n_layer-1),  # _δψ
                zeros(_n_layer-1)   # _λ_thermal
    )
);
