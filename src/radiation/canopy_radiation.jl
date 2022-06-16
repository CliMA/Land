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
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to store canopy radiation profiles

# Fields

$(TYPEDFIELDS)

"""
mutable struct BroadbandSLCanopyRadiationProfile{FT} <: AbstractCanopyRadiationProfile{FT}
    # diagnostic variables that change with time
    "Mean shaded leaf APAR (per leaf area) in μmol m⁻² s⁻¹"
    apar_shaded::FT
    "Mean sunlit leaf APAR (per leaf area) in μmol m⁻² s⁻¹"
    apar_sunlit::FT
    "Weighted extinction coefficient for diffuse radiation (ratio between projected area to true leaf area)"
    k_diffuse::FT
    "Weighted extinction coefficient for direct radiation (ratio between projected area to true leaf area)"
    k_direct::FT
    "Total shaded leaf area index"
    lai_shaded::FT
    "Total sunlit leaf area index"
    lai_sunlit::FT
    "Mean shaded leaf PAR (per leaf area) in μmol m⁻² s⁻¹"
    par_shaded::FT
    "Mean sunlit leaf PAR (per leaf area) in μmol m⁻² s⁻¹"
    par_sunlit::FT
    "Net absorbed radiation for shaded leaves `[W m⁻²]`"
    r_net_shaded::FT
    "Net absorbed radiation for sunlit leaves `[W m⁻²]`"
    r_net_sunlit::FT

    # caches to speed up calculations
    "Extinction coefficient for diffuse radiation at different leaf inclination angles"
    _k_diffuse::Vector{FT}
    "Extinction coefficient for direct radiation at different leaf inclination angles"
    _k_direct::Vector{FT}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-15: add constructor
#     2022-Jun-16: add cache values for diffuse and direct radiation
#     2022-Jun-16: add more variable to store partitions and radiations
#
#######################################################################################################################################################################################################
"""

    BroadbandSLCanopyRadiationProfile{FT}(; n_incl::Int = 9) where {FT<:AbstractFloat}

Construct a struct to store broadband canopy radiation profiles, given
- `n_incl` Number of inclination angles
"""
BroadbandSLCanopyRadiationProfile{FT}(; n_incl::Int = 9) where {FT<:AbstractFloat} = (
    return BroadbandSLCanopyRadiationProfile{FT}(
                0,                  # apar_shaded
                0,                  # apar_sunlit
                0,                  # k_diffuse
                0,                  # k_direct
                0,                  # lai_shaded
                0,                  # lai_sunlit
                0,                  # par_shaded
                0,                  # par_sunlit
                0,                  # r_net_shaded
                0,                  # r_net_sunlit
                zeros(FT,n_incl),   # _k_diffuse
                zeros(FT,n_incl)    # _k_direct
    )
);


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
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to store canopy radiation profiles

# Fields

$(TYPEDFIELDS)

"""
mutable struct HyperspectralMLCanopyRadiationProfile{FT} <: AbstractCanopyRadiationProfile{FT}
    # diagnostic variables that change with time
    "Albedo towards the viewing direction"
    albedo::Vector{FT}
    "Mean APAR for shaded leaves `[μmol m⁻² s⁻¹]`"
    apar_shaded::Vector{FT}
    "APAR for sunlit leaves `[μmol m⁻² s⁻¹]`"
    apar_sunlit::Array{FT,3}
    "Downwelling diffuse short-wave radiation at each canopy layer boundary `[mW m⁻² nm⁻¹]`"
    e_diffuse_down::Matrix{FT}
    "Upwelling diffuse short-wave radiation at each canopy layer boundary `[mW m⁻² nm⁻¹]`"
    e_diffuse_up::Matrix{FT}
    "Solar directly radiation at each canopy layer boundary `[mW m⁻² nm⁻¹]`"
    e_direct::Matrix{FT}
    "Net diffuse radiation at each canopy layer for APAR `[mW m⁻² nm⁻¹]`"
    e_net_diffuse::Matrix{FT}
    "Net direct radiation at each canopy layer for APAR `[mW m⁻² nm⁻¹]`"
    e_net_direct::Matrix{FT}
    "Total radiation towards the viewing direction `[mW m⁻² nm⁻¹]`"
    e_o::Vector{FT}
    "Sum diffuse radiation at each canopy layer for PAR `[mW m⁻² nm⁻¹]`"
    e_sum_diffuse::Matrix{FT}
    "Sum direct radiation at each canopy layer for PAR `[mW m⁻² nm⁻¹]`"
    e_sum_direct::Matrix{FT}
    "Radiation towards the viewing direction per layer (including soil) `[mW m⁻² nm⁻¹]`"
    e_v::Matrix{FT}
    "Total incoming radiation PAR `[μmol m⁻² s⁻¹]`"
    par_in::FT
    "Diffuse incoming radiation PAR `[μmol m⁻² s⁻¹]`"
    par_in_diffuse::FT
    "Direct incoming radiation PAR `[μmol m⁻² s⁻¹]`"
    par_in_direct::FT
    "Mean PAR for shaded leaves (before absorption) `[μmol m⁻² s⁻¹]`"
    par_shaded::Vector{FT}
    "PAR for sunlit leaves (before absorption) `[μmol m⁻² s⁻¹]`"
    par_sunlit::Array{FT,3}
    "Mean APAR for shaded leaves for photosynthesis `[μmol m⁻² s⁻¹]`"
    ppar_shaded::Vector{FT}
    "APAR for sunlit leaves for photosynthesis `[μmol m⁻² s⁻¹]`"
    ppar_sunlit::Array{FT,3}
    "Longwave energy flux from leaves per leaf area (one side) `[W m⁻²]`"
    r_lw::Vector{FT}
    "Downwelling longwave energy flux `[W m⁻²]`"
    r_lw_down::Vector{FT}
    "Upwelling longwave energy flux `[W m⁻²]`"
    r_lw_up::Vector{FT}
    "Net longwave energy absorption for all leaves `[W m⁻²]`"
    r_net_lw::Vector{FT}
    "Net shortwave energy absorption for all leaves `[W m⁻²]`"
    r_net_sw::Vector{FT}
    "Net shortwave energy absorption for shaded leaves `[W m⁻²]`"
    r_net_sw_shaded::Vector{FT}
    "Net shortwave energy absorption for sunlit leaves `[W m⁻²]`"
    r_net_sw_sunlit::Vector{FT}
    "Downwelling SIF for sunlit leaves at each wavelength for a layer"
    s_layer_down::Matrix{FT}
    "Upwelling SIF for sunlit leaves at each wavelength for a layer"
    s_layer_up::Matrix{FT}
    "Downwelling SIF"
    sif_down::Matrix{FT}
    "SIF at observer direction"
    sif_obs::Vector{FT}
    "SIF at observer direction from shaded APAR"
    sif_obs_shaded::Vector{FT}
    "SIF at observer direction from scattering"
    sif_obs_scatter::Vector{FT}
    "SIF at observer direction from soil reflection"
    sif_obs_ssoil::Vector{FT}
    "SIF at observer direction from sunlit APAR"
    sif_obs_sunlit::Vector{FT}
    "Upwelling SIF"
    sif_up::Matrix{FT}
    "Shaded leaf fluorescence quantum yield"
    ϕ_shaded::Vector{FT}
    "Sunlit leaf fluorescence quantum yield"
    ϕ_sunlit::Array{FT,3}

    # caches to speed up calculations
    "Mean APAR for shaded leaves per wavelength `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _apar_shaded::Vector{FT}
    "APAR for sunlit leaves per wavelength `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _apar_sunlit::Vector{FT}
    "Mean PAR for shaded leaves per wavelength (before absorption) `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _par_shaded::Vector{FT}
    "PAR for sunlit leaves per wavelength (before absorption) `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _par_sunlit::Vector{FT}
    "Mean APAR for shaded leaves for photosynthesis per wavelength `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _ppar_shaded::Vector{FT}
    "APAR for sunlit leaves for photosynthesis per wavelength `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _ppar_sunlit::Vector{FT}
    "Downwelling longwave energy flux cache `[W m⁻²]`"
    _r_emit_down::Vector{FT}
    "Upwelling longwave energy flux cache `[W m⁻²]`"
    _r_emit_up::Vector{FT}
    "Downwelling SIF for sunlit leaves at each wavelength"
    _s_emit_down::Matrix{FT}
    "Upwelling SIF for sunlit leaves at each wavelength"
    _s_emit_up::Matrix{FT}
    "Downwelling SIF for shaded leaves at each wavelength"
    _s_shaded_down::Vector{FT}
    "Upwelling SIF for shaded leaves at each wavelength"
    _s_shaded_up::Vector{FT}
    "Downwelling SIF for sunlit leaves at each wavelength"
    _s_sunlit_down::Vector{FT}
    "Upwelling SIF for sunlit leaves at each wavelength"
    _s_sunlit_up::Vector{FT}
    "Cache to compute SIF at observer direction from shaded APAR"
    _sif_obs_shaded::Matrix{FT}
    "Cache to compute SIF at observer direction from scattering"
    _sif_obs_scatter::Matrix{FT}
    "Cache to compute SIF at observer direction from sunlit APAR"
    _sif_obs_sunlit::Matrix{FT}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-09: add constructor
#     2022-Jun-09: add fields: albedo, apar_shaded, apar_sunlit, e_net_diffuse, e_net_direct, e_o, e_v, par_shaded, par_sunlit, r_net
#     2022-Jun-10: add fields: e_sum_diffuse, e_sum_direct, par_in, par_in_diffuse, par_in_direct, par_shaded, par_sunlit, _par_shaded, _par_sunlit
#     2022-Jun-10: add fields: r_net_sw, r_net_sw_shaded, r_net_sw_sunlit, r_lw, r_lw_down, r_lw_up, _r_emit_down, _r_emit_up
#     2022-Jun-10: add n_par to options and fix dimensions of the variables
#     2022-Jun-10: add n_λf for SIF
#     2022-Jun-13: add more fields for sif calculations
#     2022-Jun-15: rename to HyperspectralMLCanopyRadiationProfile
#     2022-Jun-16: fox documentation
#
#######################################################################################################################################################################################################
"""

    HyperspectralMLCanopyRadiationProfile{FT}(; n_azi::Int = 36, n_incl::Int = 9, n_layer::Int = 20, n_par::Int = 35, n_λ::Int = 114, n_λf::Int = 29) where {FT<:AbstractFloat}

Construct a struct to store hyperspectral canopy radiation profiles, given
- `n_azi` Number of azimuth angles
- `n_incl` Number of inclination angles
- `n_layer` Number of canopy layers
- `n_par` Number of PAR wavelength bins
- `n_λ` Number of wavelength bins
- `n_λf` Number of SIF wavelength bins
"""
HyperspectralMLCanopyRadiationProfile{FT}(; n_azi::Int = 36, n_incl::Int = 9, n_layer::Int = 20, n_par::Int = 35, n_λ::Int = 114, n_λf::Int = 29) where {FT<:AbstractFloat} = (
    return HyperspectralMLCanopyRadiationProfile{FT}(
                zeros(FT,n_λ),                  # albedo
                zeros(FT,n_layer),              # apar_shaded
                zeros(FT,n_incl,n_azi,n_layer), # apar_sunlit
                zeros(FT,n_λ,n_layer+1),        # e_diffuse_down
                zeros(FT,n_λ,n_layer+1),        # e_diffuse_up
                zeros(FT,n_λ,n_layer+1),        # e_direct
                zeros(FT,n_λ,n_layer),          # e_net_diffuse
                zeros(FT,n_λ,n_layer),          # e_net_direct
                zeros(FT,n_λ),                  # e_o
                zeros(FT,n_λ,n_layer),          # e_sum_diffuse
                zeros(FT,n_λ,n_layer),          # e_sum_direct
                zeros(FT,n_λ,n_layer+1),        # e_v
                0,                              # par_in
                0,                              # par_in_diffuse
                0,                              # par_in_direct
                zeros(FT,n_layer),              # par_shaded
                zeros(FT,n_incl,n_azi,n_layer), # par_sunlit
                zeros(FT,n_layer),              # ppar_shaded
                zeros(FT,n_incl,n_azi,n_layer), # ppar_sunlit
                zeros(FT,n_layer),              # r_lw
                zeros(FT,n_layer+1),            # r_lw_down
                zeros(FT,n_layer+1),            # r_lw_up
                zeros(FT,n_layer),              # r_net_lw
                zeros(FT,n_layer),              # r_net_sw
                zeros(FT,n_layer),              # r_net_sw_shaded
                zeros(FT,n_layer),              # r_net_sw_sunlit
                zeros(FT,n_λf,n_layer),         # s_layer_down
                zeros(FT,n_λf,n_layer),         # s_layer_up
                zeros(FT,n_λf,n_layer+1),       # sif_down
                zeros(FT,n_λf),                 # sif_obs
                zeros(FT,n_λf),                 # sif_obs_shaded
                zeros(FT,n_λf),                 # sif_obs_scatter
                zeros(FT,n_λf),                 # sif_obs_soil
                zeros(FT,n_λf),                 # sif_obs_sunlit
                zeros(FT,n_λf,n_layer+1),       # sif_up
                zeros(FT,n_layer),              # ϕ_shaded
                zeros(FT,n_incl,n_azi,n_layer), # ϕ_sunlit
                zeros(FT,n_par),                # _apar_shaded
                zeros(FT,n_par),                # _apar_sunlit
                zeros(FT,n_par),                # _par_shaded
                zeros(FT,n_par),                # _par_sunlit
                zeros(FT,n_par),                # _ppar_shaded
                zeros(FT,n_par),                # _ppar_sunlit
                zeros(FT,n_layer),              # _r_emit_down
                zeros(FT,n_layer+1),            # _r_emit_up
                zeros(FT,n_λf,n_layer),         # _s_emit_down
                zeros(FT,n_λf,n_layer+1),       # _s_emit_up
                zeros(FT,n_λf),                 # _s_shaded_down
                zeros(FT,n_λf),                 # _s_shaded_up
                zeros(FT,n_λf),                 # _s_sunlit_down
                zeros(FT,n_λf),                 # _s_sunlit_up
                zeros(FT,n_λf,n_layer),         # _sif_obs_shaded
                zeros(FT,n_λf,n_layer),         # _sif_obs_scatter
                zeros(FT,n_λf,n_layer)          # _sif_obs_sunlit
    )
);
