###############################################################################
#
# Cache of values to speed the canopy_fluxes!
#
###############################################################################
"""
    mutable struct CFCache{FT}

Cache to speed [`canopy_fluxes!`](@ref) by pre-allocating arrays

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct CFCache{FT}
    "absorbed energy from wave lengths"
    abs_wave   ::Array{FT,1}
    "absfs' * lidf [nAzi]"
    absfs_lidf ::Array{FT,1}
    "wave length energy [same as dWL]"
    E_all      ::Array{FT,1}
    "wave length energy [same as iPAR]"
    E_iPAR     ::Array{FT,1}
    "lPs [nLayer]"
    lPs        ::Array{FT,1}
    "kChlrel [same as iPAR]"
    kChlrel    ::Array{FT,1}
    "diffusive PAR [same as iPAR]"
    PAR_diff   ::Array{FT,1}
    "diffusive PAR for photosynthesis [same as iPAR]"
    PAR_diffCab::Array{FT,1}
    "direct PAR [same as iPAR]"
    PAR_dir    ::Array{FT,1}
    "diffusive PAR for photosynthesis [same as iPAR]"
    PAR_dirCab ::Array{FT,1}
end








###############################################################################
#
# Cache of values to speed the canopy_geometry!
#
###############################################################################
"""
    mutable struct CGCache{FT}

Cache to speed [`canopy_geometry!`](@ref) by pre-allocating arrays

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct CGCache{FT}
    # 1D arrays
    "cos_ttli .* cos(vza) dim: nIncl"
    _Co ::Array{FT,1}
    "cos_ttli .* cos(sza) dim: nIncl"
    _Cs ::Array{FT,1}
    "sin_ttli .* sin(vza) dim: nIncl"
    _So ::Array{FT,1}
    "sin_ttli .* sin(sza) dim: nIncl"
    _Ss ::Array{FT,1}

    # 2D arrays
    "maxtrix filled with 1 dim: (1, nAzi)"
    _1s ::Array{FT,2}
    "2D array to speed up _cds and _cdo dim: (nIncl, nAzi)"
    _2d ::Array{FT,2}
    "_Co * _1s .+ _So * cos_philo' dim: (nIncl, nAzi)"
    _cdo::Array{FT,2}
    "_Cs * _1s .+ _Ss * cos_ttlo' dim: (nIncl, nAzi)"
    _cds::Array{FT,2}
end








###############################################################################
#
# Cache of values to speed the short_wave!
#
###############################################################################
"""
    mutable struct SFCache{FT}

Cache to speed [`SIF_fluxes!`](@ref) by pre-allocating arrays

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct SFCache{FT}
    # 1D array
    M⁻_sun        ::Array{FT,1}
    M⁺_sun        ::Array{FT,1}
    wfEs          ::Array{FT,1}
    sfEs          ::Array{FT,1}
    sbEs          ::Array{FT,1}
    M⁺⁻           ::Array{FT,1}
    M⁺⁺           ::Array{FT,1}
    M⁻⁺           ::Array{FT,1}
    M⁻⁻           ::Array{FT,1}
    sun_dwl_iWlE  ::Array{FT,1}
    tmp_dwl_iWlE  ::Array{FT,1}
    ϕ_cosΘ_lidf   ::Array{FT,1}
    vfEplu_shade  ::Array{FT,1}
    vbEmin_shade  ::Array{FT,1}
    vfEplu_sun    ::Array{FT,1}
    vbEmin_sun    ::Array{FT,1}
    sigfEmin_shade::Array{FT,1}
    sigbEmin_shade::Array{FT,1}
    sigfEmin_sun  ::Array{FT,1}
    sigbEmin_sun  ::Array{FT,1}
    sigfEplu_shade::Array{FT,1}
    sigbEplu_shade::Array{FT,1}
    sigfEplu_sun  ::Array{FT,1}
    sigbEplu_sun  ::Array{FT,1}
    zeroB         ::Array{FT,1}
    tmp_1d_nWlF   ::Array{FT,1}
    tmp_1d_nLayer ::Array{FT,1}
    dnorm         ::Array{FT,1}

    # 2D array
    "transmission of diffusive light?"
    τ_dd                ::Array{FT,2}
    "extinction of diffuse light?"
    ρ_dd                ::Array{FT,2}
    Xdd                 ::Array{FT,2}
    Rdd                 ::Array{FT,2}
    Y                   ::Array{FT,2}
    U                   ::Array{FT,2}
    S⁻                  ::Array{FT,2}
    S⁺                  ::Array{FT,2}
    piLs                ::Array{FT,2}
    piLd                ::Array{FT,2}
    Fsmin               ::Array{FT,2}
    Fsplu               ::Array{FT,2}
    Fdmin               ::Array{FT,2}
    Fdplu               ::Array{FT,2}
    Femo                ::Array{FT,2}
    M⁺                  ::Array{FT,2}
    M⁻                  ::Array{FT,2}
    ϕ_cosΘ              ::Array{FT,2}
    F⁻                  ::Array{FT,2}
    F⁺                  ::Array{FT,2}
    net_diffuse         ::Array{FT,2}
    tmp_2d_nWlF_nLayer  ::Array{FT,2}
    tmp_2d_nWlF_nLayer_2::Array{FT,2}
end








###############################################################################
#
# Cache of values to speed the short_wave!
#
###############################################################################
"""
    mutable struct SWCache{FT}

Cache to speed [`short_wave!`](@ref) by pre-allocating arrays

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct SWCache{FT}
    # 1D arrays
    "dnorm?"
    dnorm ::Array{FT,1}
    "pi * Lo"
    piLo  ::Array{FT,1}
    "pi * Lo from canopy"
    piLoc ::Array{FT,1}
    "pi * Lo from soil"
    piLos ::Array{FT,1}

    # 2D arrays
    "pi * Lo from canopy 2D matrix"
    piLoc2::Array{FT,2}
    "extinction of diffuse light?"
    ρ_dd  ::Array{FT,2}
    "extinction of direct light?"
    ρ_sd  ::Array{FT,2}
    "transmission of diffusive light?"
    τ_dd  ::Array{FT,2}
    "transmission of direct light?"
    τ_sd  ::Array{FT,2}
end








###############################################################################
#
# Cache of values to speed up the calculations
#
###############################################################################
"""
    mutable struct RTCache{FT}

Collection of caches to speed up RT module

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct RTCache{FT}
    "[`CFCache`](@ref) type cache"
    cf_con::CFCache{FT}
    "[`CGCache`](@ref) type cache"
    cg_con::CGCache{FT}
    "[`SFCache`](@ref) type cache"
    sf_con::SFCache{FT}
    "[`SWCache`](@ref) type cache"
    sw_con::SWCache{FT}
end
