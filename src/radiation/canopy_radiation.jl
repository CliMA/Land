#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jun-15: add abstract type for canopy radiation profile
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractCanopy:
- [`BroadbandSLCanopyRadiationProfile`](@ref)
- [`HyperspectralMLCanopyRadiationProfile`](@ref)

"""
abstract type AbstractCanopyRadiationProfile{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-15: add struct for broadband radiation
#     2022-Jun-16: add cache values for diffuse and direct radiation
#     2022-Jun-16: add more variable to store partitions and radiations
#     2022-Jul-19: use kwdef for the constructor
#     2022-Jul-19: add dimension control to struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to store canopy radiation profiles

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct BroadbandSLCanopyRadiationProfile{FT<:AbstractFloat} <: AbstractCanopyRadiationProfile{FT}
    # dimensions
    "Dimension of inclination angles"
    DIM_INCL::Int = 9

    # diagnostic variables that change with time
    "Mean shaded leaf APAR (per leaf area) in μmol m⁻² s⁻¹"
    apar_shaded::FT = 0
    "Mean sunlit leaf APAR (per leaf area) in μmol m⁻² s⁻¹"
    apar_sunlit::FT = 0
    "Weighted extinction coefficient for diffuse radiation (ratio between projected area to true leaf area)"
    k_diffuse::FT = 0
    "Weighted extinction coefficient for direct radiation (ratio between projected area to true leaf area)"
    k_direct::FT = 0
    "Total shaded leaf area index"
    lai_shaded::FT = 0
    "Total sunlit leaf area index"
    lai_sunlit::FT = 0
    "Mean shaded leaf PAR (per leaf area) in μmol m⁻² s⁻¹"
    par_shaded::FT = 0
    "Mean sunlit leaf PAR (per leaf area) in μmol m⁻² s⁻¹"
    par_sunlit::FT = 0
    "Net absorbed radiation for shaded leaves `[W m⁻²]`"
    r_net_shaded::FT = 0
    "Net absorbed radiation for sunlit leaves `[W m⁻²]`"
    r_net_sunlit::FT = 0

    # caches to speed up calculations
    "Extinction coefficient for diffuse radiation at different leaf inclination angles"
    _k_diffuse::Vector{FT} = zeros(FT, DIM_INCL)
    "Extinction coefficient for direct radiation at different leaf inclination angles"
    _k_direct::Vector{FT} = zeros(FT, DIM_INCL)
end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-09: migrate CanopyRads as HyperspectralMLCanopyRadiationProfile
#     2022-Jun-09: add fields: albedo, apar_shaded, apar_sunlit, e_net_diffuse, e_net_direct, e_o, e_v, par_shaded, par_sunlit, r_net
#     2022-Jun-10: add fields: e_sum_diffuse, e_sum_direct, par_in, par_in_diffuse, par_in_direct, par_shaded, par_sunlit, _par_shaded, _par_sunlit
#     2022-Jun-10: add fields: r_net_sw, r_net_sw_shaded, r_net_sw_sunlit, r_lw, r_lw_down, r_lw_up, _r_emit_down, _r_emit_up
#     2022-Jun-10: add more fields for SIF
#     2022-Jun-13: add more fields for sif calculations
#     2022-Jun-15: rename to HyperspectralMLCanopyRadiationProfile
#     2022-Jun-27: move ppar_sunlit and ppar_shaded to Leaves2D
#     2022-Jul-19: use kwdef for the constructor
#     2022-Jul-19: add dimension control to struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to store canopy radiation profiles

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct HyperspectralMLCanopyRadiationProfile{FT<:AbstractFloat} <: AbstractCanopyRadiationProfile{FT}
    # dimensions
    "Dimension of azimuth angles"
    DIM_AZI::Int = 36
    "Dimension of inclination angles"
    DIM_INCL::Int = 9
    "Dimension of canopy layers"
    DIM_LAYER::Int = 20
    "Dimension of PAR wave length bins"
    DIM_PAR::Int = 35
    "Dimension of SIF wave length bins"
    DIM_SIF::Int = 29
    "Dimension of short wave length bins"
    DIM_WL::Int = 114

    # diagnostic variables that change with time
    "Albedo towards the viewing direction"
    albedo::Vector{FT} = zeros(FT, DIM_WL)
    "Mean APAR for shaded leaves `[μmol m⁻² s⁻¹]`"
    apar_shaded::Vector{FT} = zeros(FT, DIM_LAYER)
    "APAR for sunlit leaves `[μmol m⁻² s⁻¹]`"
    apar_sunlit::Array{FT,3} = zeros(FT, DIM_INCL, DIM_AZI, DIM_LAYER)
    "Downwelling diffuse short-wave radiation at each canopy layer boundary `[mW m⁻² nm⁻¹]`"
    e_diffuse_down::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER+1)
    "Upwelling diffuse short-wave radiation at each canopy layer boundary `[mW m⁻² nm⁻¹]`"
    e_diffuse_up::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER+1)
    "Solar directly radiation at each canopy layer boundary `[mW m⁻² nm⁻¹]`"
    e_direct::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER+1)
    "Net diffuse radiation at each canopy layer for APAR `[mW m⁻² nm⁻¹]`"
    e_net_diffuse::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Net direct radiation at each canopy layer for APAR `[mW m⁻² nm⁻¹]`"
    e_net_direct::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Total radiation towards the viewing direction `[mW m⁻² nm⁻¹]`"
    e_o::Vector{FT} = zeros(FT, DIM_WL)
    "Sum diffuse radiation at each canopy layer for PAR `[mW m⁻² nm⁻¹]`"
    e_sum_diffuse::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Sum direct radiation at each canopy layer for PAR `[mW m⁻² nm⁻¹]`"
    e_sum_direct::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Radiation towards the viewing direction per layer (including soil) `[mW m⁻² nm⁻¹]`"
    e_v::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER+1)
    "Total incoming radiation PAR `[μmol m⁻² s⁻¹]`"
    par_in::FT = 0
    "Diffuse incoming radiation PAR `[μmol m⁻² s⁻¹]`"
    par_in_diffuse::FT = 0
    "Direct incoming radiation PAR `[μmol m⁻² s⁻¹]`"
    par_in_direct::FT = 0
    "Mean PAR for shaded leaves (before absorption) `[μmol m⁻² s⁻¹]`"
    par_shaded::Vector{FT} = zeros(FT, DIM_LAYER)
    "PAR for sunlit leaves (before absorption) `[μmol m⁻² s⁻¹]`"
    par_sunlit::Array{FT,3} = zeros(FT, DIM_INCL, DIM_AZI, DIM_LAYER)
    "Longwave energy flux from leaves per leaf area (one side) `[W m⁻²]`"
    r_lw::Vector{FT} = zeros(FT, DIM_LAYER)
    "Downwelling longwave energy flux `[W m⁻²]`"
    r_lw_down::Vector{FT} = zeros(FT, DIM_LAYER+1)
    "Upwelling longwave energy flux `[W m⁻²]`"
    r_lw_up::Vector{FT} = zeros(FT, DIM_LAYER+1)
    "Net longwave energy absorption for all leaves `[W m⁻²]`"
    r_net_lw::Vector{FT} = zeros(FT, DIM_LAYER)
    "Net shortwave energy absorption for all leaves `[W m⁻²]`"
    r_net_sw::Vector{FT} = zeros(FT, DIM_LAYER)
    "Net shortwave energy absorption for shaded leaves `[W m⁻²]`"
    r_net_sw_shaded::Vector{FT} = zeros(FT, DIM_LAYER)
    "Net shortwave energy absorption for sunlit leaves `[W m⁻²]`"
    r_net_sw_sunlit::Vector{FT} = zeros(FT, DIM_LAYER)
    "Downwelling SIF for sunlit leaves at each wavelength for a layer"
    s_layer_down::Matrix{FT} = zeros(FT, DIM_SIF, DIM_LAYER)
    "Upwelling SIF for sunlit leaves at each wavelength for a layer"
    s_layer_up::Matrix{FT} = zeros(FT, DIM_SIF, DIM_LAYER)
    "Downwelling SIF"
    sif_down::Matrix{FT} = zeros(FT, DIM_SIF, DIM_LAYER+1)
    "SIF at observer direction"
    sif_obs::Vector{FT} = zeros(FT, DIM_SIF)
    "SIF at observer direction from shaded APAR"
    sif_obs_shaded::Vector{FT} = zeros(FT, DIM_SIF)
    "SIF at observer direction from scattering"
    sif_obs_scatter::Vector{FT} = zeros(FT, DIM_SIF)
    "SIF at observer direction from soil reflection"
    sif_obs_ssoil::Vector{FT} = zeros(FT, DIM_SIF)
    "SIF at observer direction from sunlit APAR"
    sif_obs_sunlit::Vector{FT} = zeros(FT, DIM_SIF)
    "Upwelling SIF"
    sif_up::Matrix{FT} = zeros(FT, DIM_SIF, DIM_LAYER+1)

    # caches to speed up calculations
    "Mean APAR for shaded leaves per wavelength `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _apar_shaded::Vector{FT} = zeros(FT, DIM_PAR)
    "APAR for sunlit leaves per wavelength `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _apar_sunlit::Vector{FT} = zeros(FT, DIM_PAR)
    "Mean PAR for shaded leaves per wavelength (before absorption) `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _par_shaded::Vector{FT} = zeros(FT, DIM_PAR)
    "PAR for sunlit leaves per wavelength (before absorption) `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _par_sunlit::Vector{FT} = zeros(FT, DIM_PAR)
    "Mean APAR for shaded leaves for photosynthesis per wavelength `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _ppar_shaded::Vector{FT} = zeros(FT, DIM_PAR)
    "APAR for sunlit leaves for photosynthesis per wavelength `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _ppar_sunlit::Vector{FT} = zeros(FT, DIM_PAR)
    "Downwelling longwave energy flux cache `[W m⁻²]`"
    _r_emit_down::Vector{FT} = zeros(FT, DIM_LAYER)
    "Upwelling longwave energy flux cache `[W m⁻²]`"
    _r_emit_up::Vector{FT} = zeros(FT, DIM_LAYER+1)
    "Downwelling SIF for sunlit leaves at each wavelength"
    _s_emit_down::Matrix{FT} = zeros(FT, DIM_SIF, DIM_LAYER)
    "Upwelling SIF for sunlit leaves at each wavelength"
    _s_emit_up::Matrix{FT} = zeros(FT, DIM_SIF, DIM_LAYER+1)
    "Downwelling SIF for shaded leaves at each wavelength"
    _s_shaded_down::Vector{FT} = zeros(FT, DIM_SIF)
    "Upwelling SIF for shaded leaves at each wavelength"
    _s_shaded_up::Vector{FT} = zeros(FT, DIM_SIF)
    "Downwelling SIF for sunlit leaves at each wavelength"
    _s_sunlit_down::Vector{FT} = zeros(FT, DIM_SIF)
    "Upwelling SIF for sunlit leaves at each wavelength"
    _s_sunlit_up::Vector{FT} = zeros(FT, DIM_SIF)
    "Cache to compute SIF at observer direction from shaded APAR"
    _sif_obs_shaded::Matrix{FT} = zeros(FT, DIM_SIF, DIM_LAYER)
    "Cache to compute SIF at observer direction from scattering"
    _sif_obs_scatter::Matrix{FT} = zeros(FT, DIM_SIF, DIM_LAYER)
    "Cache to compute SIF at observer direction from sunlit APAR"
    _sif_obs_sunlit::Matrix{FT} = zeros(FT, DIM_SIF, DIM_LAYER)
end
