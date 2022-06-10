#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-09: migrate CanopyRads as CanopyRadiationProfile
#     2022-Jun-09: add fields: albedo, apar_shaded, apar_sunlit, e_net_diffuse, e_net_direct, e_o, e_v, par_shaded, par_sunlit, r_net
#     2022-Jun-10: add fields: e_sum_diffuse, e_sum_direct, par_in, par_in_diffuse, par_in_direct, par_shaded, par_sunlit, _par_shaded, _par_sunlit
#     2022-Jun-10: add fields: r_net_sw, r_net_sw_shaded, r_net_sw_sunlit
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
    "Net energy absorption for all leaves `[W m⁻²]`"
    r_net_sw::Vector{FT}
    "Net energy absorption for shaded leaves `[W m⁻²]`"
    r_net_sw_shaded::Vector{FT}
    "Net energy absorption for sunlit leaves `[W m⁻²]`"
    r_net_sw_sunlit::Vector{FT}
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
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-09: add constructor
#     2022-Jun-09: add fields: albedo, apar_shaded, apar_sunlit, e_net_diffuse, e_net_direct, e_o, e_v, par_shaded, par_sunlit, r_net
#     2022-Jun-10: add fields: e_sum_diffuse, e_sum_direct, par_in, par_in_diffuse, par_in_direct, par_shaded, par_sunlit, _par_shaded, _par_sunlit
#     2022-Jun-10: add fields: r_net_sw, r_net_sw_shaded, r_net_sw_sunlit
#     2022-Jun-10: add n_par to options and fix dimensions of the variables
#
#######################################################################################################################################################################################################
"""

    CanopyRadiationProfile{FT}(; n_azi::Int = 36, n_incl::Int = 9, n_layer::Int = 20, n_par::Int = 35, n_λ::Int = 114) where {FT<:AbstractFloat}

Construct a struct to store canopy radiation profiles, given
- `n_azi` Number of azimuth angles
- `n_incl` Number of inclination angles
- `n_layer` Number of canopy layers
- `n_par` Number of PAR wavelength bins
- `n_λ` Number of wavelength bins
"""
CanopyRadiationProfile{FT}(; n_azi::Int = 36, n_incl::Int = 9, n_layer::Int = 20, n_par::Int = 35, n_λ::Int = 114) where {FT<:AbstractFloat} = (
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
                zeros(FT,n_layer),              # r_net_sw
                zeros(FT,n_layer),              # r_net_sw_shaded
                zeros(FT,n_layer),              # r_net_sw_sunlit
                zeros(FT,n_par),                # _apar_shaded
                zeros(FT,n_par),                # _apar_sunlit
                zeros(FT,n_par),                # _par_shaded
                zeros(FT,n_par),                # _par_sunlit
                zeros(FT,n_par),                # _ppar_shaded
                zeros(FT,n_par)                 # _ppar_sunlit
    )
);
