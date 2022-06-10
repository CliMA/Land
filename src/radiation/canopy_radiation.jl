#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-09: migrate CanopyRads as CanopyRadiationProfile
#     2022-Jun-09: add fields: albedo, apar_shaded, apar_sunlit, e_net_diffuse, e_net_direct, e_o, e_v, par_shaded, par_sunlit, r_net
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to store canopy radiation profiles

# Fields

$(TYPEDFIELDS)

"""
mutable struct CanopyRadiationProfile{FT<:AbstractFloat}
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
    "Net diffuse radiation at each canopy layer `[mW m⁻² nm⁻¹]`"
    e_net_diffuse::Matrix{FT}
    "Net direct radiation at each canopy layer `[mW m⁻² nm⁻¹]`"
    e_net_direct::Matrix{FT}
    "Total radiation towards the viewing direction `[mW m⁻² nm⁻¹]`"
    e_o::Vector{FT}
    "Radiation towards the viewing direction per layer (including soil) `[mW m⁻² nm⁻¹]`"
    e_v::Matrix{FT}
    "Mean APAR for shaded leaves for photosynthesis `[μmol m⁻² s⁻¹]`"
    ppar_shaded::Vector{FT}
    "APAR for sunlit leaves for photosynthesis `[μmol m⁻² s⁻¹]`"
    ppar_sunlit::Array{FT,3}
    "Net energy absorption for shortwave `[W m⁻²]`"
    r_net::FT
    "Mean APAR for shaded leaves per wavelength `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _apar_shaded::Vector{FT}
    "APAR for sunlit leaves per wavelength `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _apar_sunlit::Vector{FT}
    "Mean APAR for shaded leaves for photosynthesis per wavelength `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _ppar_shaded::Vector{FT}
    "APAR for sunlit leaves for photosynthesis per wavelength `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _ppar_sunlit::Vector{FT}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-09: add constructor
#     2022-Jun-09: add fields: albedo, apar_shaded, apar_sunlit, e_net_diffuse, e_net_direct, e_o, e_v, par_shaded, par_sunlit, r_net
#
#######################################################################################################################################################################################################
"""

    CanopyRadiationProfile{FT}(; n_azi::Int = 36, n_incl::Int = 9, n_layer::Int = 20, n_λ::Int = 114) where {FT<:AbstractFloat}

Construct a struct to store canopy radiation profiles, given
    - `n_azi` Number of azimuth angles
    - `n_incl` Number of inclination angles
    - `n_layer` Number of canopy layers
    - `n_λ` Number of wavelength bins
"""
CanopyRadiationProfile{FT}(; n_azi::Int = 36, n_incl::Int = 9, n_layer::Int = 20, n_λ::Int = 114) where {FT<:AbstractFloat} = (
    return CanopyRadiationProfile{FT}(
                zeros(FT,n_λ),                  # albedo
                zeros(FT,n_layer),              # apar_shaded
                zeros(FT,n_incl,n_azi,n_layer), # apar_sunlit
                zeros(FT,n_λ,n_layer+1),        # e_diffuse_down
                zeros(FT,n_λ,n_layer+1),        # e_diffuse_up
                zeros(FT,n_λ,n_layer+1),        # e_direct
                zeros(FT,n_λ,n_layer),          # e_net_diffuse
                zeros(FT,n_λ,n_layer),          # e_net_direct
                zeros(FT,n_λ),                  # e_o
                zeros(FT,n_λ,n_layer+1),        # e_v
                zeros(FT,n_layer),              # ppar_shaded
                zeros(FT,n_incl,n_azi,n_layer), # ppar_sunlit
                0,                              # r_net
                zeros(FT,n_layer),              # _apar_shaded
                zeros(FT,n_layer),              # _apar_sunlit
                zeros(FT,n_layer),              # _ppar_shaded
                zeros(FT,n_layer)               # _ppar_sunlit
    )
);
