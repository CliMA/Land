###############################################################################
#
# Container of values to speed up the calculations
# Need to move many variables to other structures
#
###############################################################################
"""
    mutable struct RTContainer{FT}

A struct for canopy radiation information

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct RTContainer{FT}
    # these names looks very meaningless...
    "cos leaf azimuth angles"
    cos_ttlo ::Array{FT,1}
    "cos leaf azimuth angles - psi"
    cos_philo::Array{FT,1}
    "cos of normal of upperside of leaf"
    cos_ttli ::Array{FT,1}
    "sine of normal of upperside of leaf"
    sin_ttli ::Array{FT,1}
    "voscatt result"
    vol_scatt::Array{FT,1}

    # containers for more angles to use in canopy_geometry!
    "_Cs = cos_ttli .* cos(tts); # [nli]"
    _Cs::Array{FT,1}
    "_Ss = sin_ttli .* sin(tts); # [nli]"
    _Ss::Array{FT,1}
    "_Co = cos_ttli .* cos(tto); # [nli]"
    _Co::Array{FT,1}
    "_So = sin_ttli .* sin(tto); # [nli]"
    _So::Array{FT,1}
    "_1s = ones(FT,1,length(lazitab));"
    _1s::Array{FT,2}
    "cds = _Cs * _1s .+ _Ss * cos_ttlo';  # [nli, nlazi]"
    cds::Array{FT,2}
    "cdo = _Co * _1s .+ _So * cos_philo'; # [nli, nlazi]"
    cdo::Array{FT,2}
    "2D array to speed up; # [nli, nlazi]"
    _2d::Array{FT,2}

    # containers to speed short_wave!
    "transmission of diffusive light?"
    τ_dd  ::Array{FT,2}
    "transmission of direct light?"
    τ_sd  ::Array{FT,2}
    "extinction of diffuse light?"
    ρ_dd  ::Array{FT,2}
    "extinction of direct light?"
    ρ_sd  ::Array{FT,2}
    "dnorm?"
    dnorm ::Array{FT,1}
    "pi * Lo"
    piLo  ::Array{FT,1}
    "pi * Lo from canopy 2D matrix"
    piLoc2::Array{FT,2}
    "pi * Lo from canopy"
    piLoc ::Array{FT,1}
    "pi * Lo from soil"
    piLos ::Array{FT,1}

    # containers to speed up canopy_fluxes!
    "absorbed energy from wave lengths"
    abs_wave   ::Array{FT,1}
    "absfs' * lidf [nAzi]"
    absfs_lidf ::Array{FT,1}
    "lPs [nLayer]"
    lPs        ::Array{FT,1}
    "kChlrel [same as iPAR]"
    kChlrel    ::Array{FT,1}
    "wave lengths of iPAR, set once ans use forever [same as iPAR]"
    λ_iPAR     ::Array{FT,1}
    "absorbed iPAR [same as iPAR]"
    dλ_iPAR    ::Array{FT,1}
    "absorbed iPAR [same as iWlE]"
    dλ_iWlE    ::Array{FT,1}
    "wave length energy [same as iPAR]"
    E_iPAR     ::Array{FT,1}
    "wave length energy [same as dwl]"
    E_all      ::Array{FT,1}
    "diffusive PAR [same as iPAR]"
    PAR_diff   ::Array{FT,1}
    "direct PAR [same as iPAR]"
    PAR_dir    ::Array{FT,1}
    "diffusive PAR for photosynthesis [same as iPAR]"
    PAR_diffCab::Array{FT,1}
    "diffusive PAR for photosynthesis [same as iPAR]"
    PAR_dirCab ::Array{FT,1}

    # containers to speed up sif_fluxes!
    "soil albedo for SIF related wave length"
    soil_sif_albedo::Array{FT,1}
    "transmission of diffusive light?"
    τ_dd_sif       ::Array{FT,2}
    "extinction of diffuse light?"
    ρ_dd_sif       ::Array{FT,2}
    "temporary variables for sif_fluxes!"
    S⁻             ::Array{FT,2}
    "temporary variables for sif_fluxes!"
    S⁺             ::Array{FT,2}
    "temporary variables for sif_fluxes!"
    piLs           ::Array{FT,2}
    "temporary variables for sif_fluxes!"
    piLd           ::Array{FT,2}
    "temporary variables for sif_fluxes!"
    Fsmin          ::Array{FT,2}
    "temporary variables for sif_fluxes!"
    Fsplu          ::Array{FT,2}
    "temporary variables for sif_fluxes!"
    Fdmin          ::Array{FT,2}
    "temporary variables for sif_fluxes!"
    Fdplu          ::Array{FT,2}
    "temporary variables for sif_fluxes!"
    Femo           ::Array{FT,2}
    "temporary variables for sif_fluxes!"
    M⁺             ::Array{FT,2}
    "temporary variables for sif_fluxes!"
    M⁻             ::Array{FT,2}
    "temporary variables for sif_fluxes!"
    M⁻_sun         ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    M⁺_sun         ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    wfEs           ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    sfEs           ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    sbEs           ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    M⁺⁻            ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    M⁺⁺            ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    M⁻⁺            ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    M⁻⁻            ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    sun_dwl_iWlE   ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    tmp_dwl_iWlE   ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    ϕ_cosΘ         ::Array{FT,2}
    "temporary variables for sif_fluxes!"
    ϕ_cosΘ_lidf    ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    vfEplu_shade   ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    vbEmin_shade   ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    vfEplu_sun     ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    vbEmin_sun     ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    sigfEmin_shade ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    sigbEmin_shade ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    sigfEmin_sun   ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    sigbEmin_sun   ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    sigfEplu_shade ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    sigbEplu_shade ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    sigfEplu_sun   ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    sigbEplu_sun   ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    zeroB          ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    F⁻             ::Array{FT,2}
    "temporary variables for sif_fluxes!"
    F⁺             ::Array{FT,2}
    "temporary variables for sif_fluxes!"
    net_diffuse    ::Array{FT,2}
    "temporary variables for sif_fluxes!"
    tmp_1d_nWlF    ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    tmp_1d_nLayer  ::Array{FT,1}
    "temporary variables for sif_fluxes!"
    tmp_2d_nWlF_nLayer::Array{FT,2}
    "temporary variables for sif_fluxes!"
    tmp_2d_nWlF_nLayer_2::Array{FT,2}
end
