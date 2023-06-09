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
    abs_wave::Vector{FT}
    "absfs' * lidf [nAzi]"
    absfs_lidf::Vector{FT}
    "wave length energy [same as dWL]"
    E_all::Vector{FT}
    "wave length energy [same as iPAR]"
    E_iPAR::Vector{FT}
    "lPs [nLayer]"
    lPs::Vector{FT}
    "kChlrel [same as iPAR]"
    kChlrel::Vector{FT}
    "diffusive PAR [same as iPAR]"
    PAR_diff::Vector{FT}
    "diffusive PAR for photosynthesis [same as iPAR]"
    PAR_diffCab::Vector{FT}
    "direct PAR [same as iPAR]"
    PAR_dir::Vector{FT}
    "diffusive PAR for photosynthesis [same as iPAR]"
    PAR_dirCab::Vector{FT}
end

CFCache{FT}(rt_dim::RTDimensions) where {FT} = (
    (; nAzi, nLayer, nPAR, nWL) = rt_dim;

    return CFCache{FT}(
                abs_wave    = zeros(FT, nWL),
                absfs_lidf  = zeros(FT, nAzi),
                E_all       = zeros(FT, nWL),
                E_iPAR      = zeros(FT, nPAR),
                lPs         = zeros(FT, nLayer),
                kChlrel     = zeros(FT, nPAR),
                PAR_diff    = zeros(FT, nPAR),
                PAR_diffCab = zeros(FT, nPAR),
                PAR_dir     = zeros(FT, nPAR),
                PAR_dirCab  = zeros(FT, nPAR)
    )
);


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
    _Co::Vector{FT}
    "cos_ttli .* cos(sza) dim: nIncl"
    _Cs::Vector{FT}
    "sin_ttli .* sin(vza) dim: nIncl"
    _So::Vector{FT}
    "sin_ttli .* sin(sza) dim: nIncl"
    _Ss::Vector{FT}

    # 2D arrays
    "maxtrix filled with 1 dim: (1, nAzi)"
    _1s::Matrix{FT}
    "2D array to speed up _cds and _cdo dim: (nIncl, nAzi)"
    _2d::Matrix{FT}
    "_Co * _1s .+ _So * cos_philo' dim: (nIncl, nAzi)"
    _cdo::Matrix{FT}
    "_Cs * _1s .+ _Ss * cos_ttlo' dim: (nIncl, nAzi)"
    _cds::Matrix{FT}
end

CGCache{FT}(rt_dim::RTDimensions) where {FT} = (
    (; nAzi, nIncl) = rt_dim;

    return CGCache{FT}(
                _Co  = zeros(FT, nIncl),
                _Cs  = zeros(FT, nIncl),
                _So  = zeros(FT, nIncl),
                _Ss  = zeros(FT, nIncl),
                _1s  = ones(FT, (1, nAzi)),
                _2d  = zeros(FT, (nIncl, nAzi)),
                _cdo = zeros(FT, (nIncl, nAzi)),
                _cds = zeros(FT, (nIncl, nAzi))
    )
);


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
    M⁻_sun::Vector{FT}
    M⁺_sun::Vector{FT}
    wfEs::Vector{FT}
    sfEs::Vector{FT}
    sbEs::Vector{FT}
    M⁺⁻::Vector{FT}
    M⁺⁺::Vector{FT}
    M⁻⁺::Vector{FT}
    M⁻⁻::Vector{FT}
    sun_dwl_iWlE::Vector{FT}
    tmp_dwl_iWlE::Vector{FT}
    ϕ_cosΘ_lidf::Vector{FT}
    vfEplu_shade::Vector{FT}
    vbEmin_shade::Vector{FT}
    vfEplu_sun::Vector{FT}
    vbEmin_sun::Vector{FT}
    sigfEmin_shade::Vector{FT}
    sigbEmin_shade::Vector{FT}
    sigfEmin_sun::Vector{FT}
    sigbEmin_sun::Vector{FT}
    sigfEplu_shade::Vector{FT}
    sigbEplu_shade::Vector{FT}
    sigfEplu_sun::Vector{FT}
    sigbEplu_sun::Vector{FT}
    zeroB::Vector{FT}
    tmp_1d_nWlF::Vector{FT}
    tmp_1d_nLayer::Vector{FT}
    dnorm::Vector{FT}

    # 2D array
    "transmission of diffusive light?"
    τ_dd::Matrix{FT}
    "extinction of diffuse light?"
    ρ_dd::Matrix{FT}
    Xdd::Matrix{FT}
    Rdd::Matrix{FT}
    Y::Matrix{FT}
    U::Matrix{FT}
    S⁻::Matrix{FT}
    S⁺::Matrix{FT}
    piLs::Matrix{FT}
    piLd::Matrix{FT}
    Fsmin::Matrix{FT}
    Fsplu::Matrix{FT}
    Fdmin::Matrix{FT}
    Fdplu::Matrix{FT}
    Femo::Matrix{FT}
    M⁺::Matrix{FT}
    M⁻::Matrix{FT}
    ϕ_cosΘ::Matrix{FT}
    F⁻::Matrix{FT}
    F⁺::Matrix{FT}
    net_diffuse::Matrix{FT}
    tmp_2d_nWlF_nLayer::Matrix{FT}
    tmp_2d_nWlF_nLayer_2::Matrix{FT}
end

SFCache{FT}(rt_dim::RTDimensions) where {FT} = (
    (; nAzi, nIncl, nLayer, nLevel, nWLE, nWLF) = rt_dim;

    return SFCache{FT}(
                τ_dd                 = zeros(FT, (nWLF,nLayer)),
                Xdd                  = zeros(FT, (nWLF,nLayer)),
                Rdd                  = zeros(FT, (nWLF,nLevel)),
                Y                    = zeros(FT, (nWLF,nLayer)),
                U                    = zeros(FT, (nWLF,nLevel)),
                dnorm                = zeros(FT, nWLF),
                ρ_dd                 = zeros(FT, (nWLF,nLayer)),
                S⁻                   = zeros(FT, (nWLF,nLayer)),
                S⁺                   = zeros(FT, (nWLF,nLayer)),
                piLs                 = zeros(FT, (nWLF,nLayer)),
                piLd                 = zeros(FT, (nWLF,nLayer)),
                Fsmin                = zeros(FT, (nWLF,nLayer)),
                Fsplu                = zeros(FT, (nWLF,nLayer)),
                Fdmin                = zeros(FT, (nWLF,nLayer)),
                Fdplu                = zeros(FT, (nWLF,nLayer)),
                Femo                 = zeros(FT, (nWLF,nLayer)),
                M⁺                   = zeros(FT, (nWLF,nWLE)),
                M⁻                   = zeros(FT, (nWLF,nWLE)),
                M⁻_sun               = zeros(FT, nWLF),
                M⁺_sun               = zeros(FT, nWLF),
                wfEs                 = zeros(FT, nWLF),
                sfEs                 = zeros(FT, nWLF),
                sbEs                 = zeros(FT, nWLF),
                M⁺⁻                  = zeros(FT, nWLF),
                M⁺⁺                  = zeros(FT, nWLF),
                M⁻⁺                  = zeros(FT, nWLF),
                M⁻⁻                  = zeros(FT, nWLF),
                sun_dwl_iWlE         = zeros(FT, nWLE),
                tmp_dwl_iWlE         = zeros(FT, nWLE),
                ϕ_cosΘ               = zeros(FT, (nIncl,nAzi)),
                ϕ_cosΘ_lidf          = zeros(FT, nAzi),
                vfEplu_shade         = zeros(FT, nWLF),
                vbEmin_shade         = zeros(FT, nWLF),
                vfEplu_sun           = zeros(FT, nWLF),
                vbEmin_sun           = zeros(FT, nWLF),
                sigfEmin_shade       = zeros(FT, nWLF),
                sigbEmin_shade       = zeros(FT, nWLF),
                sigfEmin_sun         = zeros(FT, nWLF),
                sigbEmin_sun         = zeros(FT, nWLF),
                sigfEplu_shade       = zeros(FT, nWLF),
                sigbEplu_shade       = zeros(FT, nWLF),
                sigfEplu_sun         = zeros(FT, nWLF),
                sigbEplu_sun         = zeros(FT, nWLF),
                zeroB                = zeros(FT, nWLF),
                F⁻                   = zeros(FT, (nWLF,nLevel)),
                F⁺                   = zeros(FT, (nWLF,nLevel)),
                net_diffuse          = zeros(FT, (nWLF,nLayer)),
                tmp_1d_nWlF          = zeros(FT, nWLF),
                tmp_1d_nLayer        = zeros(FT, nLayer),
                tmp_2d_nWlF_nLayer   = zeros(FT, (nWLF,nLayer)),
                tmp_2d_nWlF_nLayer_2 = zeros(FT, (nWLF,nLayer))
    )
);


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
    dnorm::Vector{FT}
    "pi * Lo"
    piLo::Vector{FT}
    "pi * Lo from canopy"
    piLoc::Vector{FT}
    "pi * Lo from soil"
    piLos::Vector{FT}

    # 2D arrays
    "pi * Lo from canopy 2D matrix"
    piLoc2::Matrix{FT}
    "extinction of diffuse light?"
    ρ_dd::Matrix{FT}
    "extinction of direct light?"
    ρ_sd::Matrix{FT}
    "transmission of diffusive light?"
    τ_dd::Matrix{FT}
    "transmission of direct light?"
    τ_sd::Matrix{FT}
end

SWCache{FT}(rt_dim::RTDimensions) where {FT} = (
    (; nLayer, nWL) = rt_dim;

    return SWCache{FT}(
                dnorm  = zeros(FT, nWL),
                piLo   = zeros(FT, nWL),
                piLoc  = zeros(FT, nWL),
                piLos  = zeros(FT, nWL),
                piLoc2 = zeros(FT, (nWL,nLayer)),
                ρ_dd   = zeros(FT, (nWL,nLayer)),
                ρ_sd   = zeros(FT, (nWL,nLayer)),
                τ_dd   = zeros(FT, (nWL,nLayer)),
                τ_sd   = zeros(FT, (nWL,nLayer))
    )
);


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

function RTCache{FT}(rt_dim::RTDimensions) where {FT}
    return RTCache{FT}(
                cf_con = CFCache{FT}(rt_dim),
                cg_con = CGCache{FT}(rt_dim),
                sf_con = SFCache{FT}(rt_dim),
                sw_con = SWCache{FT}(rt_dim)
    )
end
