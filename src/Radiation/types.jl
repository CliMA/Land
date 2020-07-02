###############################################################################
#
# Canopy4RT
#
###############################################################################
"""
    struct Canopy4RT

A canopy struct for the radiation transfer module

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct Canopy4RT{FT<:AbstractFloat}
    "Number of canopy layers"
    nlayers   ::Int = FT(20.0)
    "Leaf Area Index"
    LAI       ::FT  = FT(3.0 )
    "Clumping factor"
    Ω         ::FT  = FT(1.0 )
    "Structure factor a"
    clump_a   ::FT  = FT(1.0 )
    "Structure factor b"
    clump_b   ::FT  = FT(0.0 )
    "Leaf width"
    leaf_width::FT  = FT(0.1 )
    "Vegetation height"
    hc        ::FT  = FT(2.0 )
    "Leaf Inclination"
    LIDFa     ::FT  = FT(0.0 )
    "Variation in leaf inclination"
    LIDFb     ::FT  = FT(0.0 )
    "HotSpot parameter (still need to check!)"
    hot       ::FT  = FT(0.05)

    # tree/canopy/leaf traits
    "Canopy height `[m]`"
    height::FT = FT(20.0  )
    "Canopy roughness `[m]`"
    z0m   ::FT = FT(1.0   )
    "Tree roughtnes `[m]`"
    z0h   ::FT = FT(-999.0)
    "Canopy displacement height `[m]`"
    d     ::FT = FT(-999.0)
    "m/sqrt(s) turbulent transfer coefficient"
    Cd    ::FT = FT(0.01  )


    # Some more derived parameters:
    litab    ::Array{FT,1} = collect(FT,5:10:85)
    litab_bnd::Array{FT,2} = [collect(0:10:80) collect(FT,10:10:90)]
    lazitab  ::Array{FT,1} = collect(FT,5:10:355)
    # This is changed afterwards, ignore here.
    lidf     ::Array{FT,1} = litab .* 0
    xl       ::Array{FT,1} = collect(0.0:-1.0/nlayers:-1.0)
    dx       ::FT          = 1.0/nlayers
end








###############################################################################
#
# Canopy optical
#
###############################################################################
"""
    AbstractCanopyOpti

A abstract type for canopy optical properties
"""
abstract type AbstractCanopyOpti end




"""
    struct CanopyOptiArray{FT, nWL, nLayer, nLevel, nAzi, nIncl}

A struct for canopy optical properties using Array

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct CanopyOptiArray{FT, nWL, nLayer, nLevel, nAzi, nIncl} <: AbstractCanopyOpti
    # Single value
    "Solar -> Diffuse backscatter weight"
    sdb::FT = FT(0.0)
    "Solar -> Diffuse forward scatter weight"
    sdf::FT = FT(0.0)
    "Diffuse -> Directional backscatter weight"
    dob::FT = FT(0.0)
    "Diffuse -> Directional forward scatter weight"
    dof::FT = FT(0.0)
    "Diffuse -> Diffuse backscatter weight"
    ddb::FT = FT(0.0)
    "Diffuse -> Diffuse forward scatter weight"
    ddf::FT = FT(0.0)
    "Solar beam extinction coefficient weight"
    ks ::FT = FT(0.0)
    "Outgoing beam extinction coefficient weight"
    ko ::FT = FT(0.0)
    # TODO What is this?
    "?"
    bf ::FT = FT(0.0)
    "Weight of specular2directional backscatter coefficient"
    sob::FT = FT(0.0)
    "Weight of specular2directional forward coefficient"
    sof::FT = FT(0.0)

    # dimension of nLevel
    "Probability of directly viewing a leaf in solar direction"
    Ps ::Array{FT,1} = zeros(FT, nLevel)
    "Probability of directly viewing a leaf in viewing direction"
    Po ::Array{FT,1} = zeros(FT, nLevel)
    "Bi-directional probability of directly viewing a leaf (solar->canopy->viewing)"
    Pso::Array{FT,1} = zeros(FT, nLevel)

    # dimension of nIncl * nAzi
    "conversion factor fs to compute irradiance on inclined leaf"
    fs     ::Array{FT,2} = zeros(FT, (nIncl, nAzi))
    "abs(fs)"
    absfs  ::Array{FT,2} = zeros(FT, (nIncl, nAzi))
    "abs(fs*fo)"
    absfsfo::Array{FT,2} = zeros(FT, (nIncl, nAzi))
    "fs*fo"
    fsfo   ::Array{FT,2} = zeros(FT, (nIncl, nAzi))
    "conversion factor fo for angle towards observer (not sun like fs)"
    fo     ::Array{FT,2} = zeros(FT, (nIncl, nAzi))
    "Cosine of leaf azimuths"
    cosΘ_l ::Array{FT,2} = zeros(FT, (nIncl, nAzi))
    "cos of leaf azimuth sqared"
    cos2Θ_l::Array{FT,2} = zeros(FT, (nIncl, nAzi))

    # The following also depend on leaf reflectance and transmission.
    # Might go into a separate strcuture so that we can have it separately for thermal, SW and SIF?
    # dimension of nWL * nLayer
    "diffuse     backscatter scattering coefficient for diffuse  incidence"
    sigb::Array{FT,2} = zeros(FT, (nWL,nLayer))
    "diffuse     forward     scattering coefficient for diffuse  incidence"
    sigf::Array{FT,2} = zeros(FT, (nWL,nLayer))
    "diffuse     backscatter scattering coefficient for specular incidence"
    sb  ::Array{FT,2} = zeros(FT, (nWL,nLayer))
    "diffuse     forward     scattering coefficient for specular incidence"
    sf  ::Array{FT,2} = zeros(FT, (nWL,nLayer))
    "directional backscatter scattering coefficient for diffuse  incidence"
    vb  ::Array{FT,2} = zeros(FT, (nWL,nLayer))
    "directional forward     scattering coefficient for diffuse  incidence"
    vf  ::Array{FT,2} = zeros(FT, (nWL,nLayer))
    "bidirectional scattering coefficent (directional-directional)"
    w   ::Array{FT,2} = zeros(FT, (nWL,nLayer))
    "attenuation"
    a   ::Array{FT,2} = zeros(FT, (nWL,nLayer))
    "Effective layer transmittance (direct->diffuse)"
    Xsd ::Array{FT,2} = zeros(FT, (nWL,nLayer))
    "Effective layer transmittance (diffuse->diffuse)"
    Xdd ::Array{FT,2} = zeros(FT, (nWL,nLayer))

    # dimension of nWL * nLevel
    "Effective layer reflectance (direct->diffuse)"
    R_sd::Array{FT,2} = zeros(FT, (nWL,nLevel))
    "Effective layer reflectance (diffuse->diffuse)"
    R_dd::Array{FT,2} = zeros(FT, (nWL,nLevel))
    "Solar direct radiation per layer)"
    Es_::Array{FT,2}  = zeros(FT, (nWL,nLevel))
end




"""
    struct CanopyOptiMArray{FT, nWL, nLayer, nLevel, nAzi, nIncl, nLayer_nWL, nLev_nWL, nIncl_nAzi}

A struct for canopy optical properties using MArray

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct CanopyOptiMArray{FT, nWL, nLayer, nLevel, nAzi, nIncl, nLayer_nWL, nLev_nWL, nIncl_nAzi} <: AbstractCanopyOpti
    # Single value
    "Solar -> Diffuse backscatter weight"
    sdb::FT = FT(0.0)
    "Solar -> Diffuse forward scatter weight"
    sdf::FT = FT(0.0)
    "Diffuse -> Directional backscatter weight"
    dob::FT = FT(0.0)
    "Diffuse -> Directional forward scatter weight"
    dof::FT = FT(0.0)
    "Diffuse -> Diffuse backscatter weight"
    ddb::FT = FT(0.0)
    "Diffuse -> Diffuse forward scatter weight"
    ddf::FT = FT(0.0)
    "Solar beam extinction coefficient weight"
    ks ::FT = FT(0.0)
    "Outgoing beam extinction coefficient weight"
    ko ::FT = FT(0.0)
    # TODO What is this?
    "?"
    bf ::FT = FT(0.0)
    "Weight of specular2directional backscatter coefficient"
    sob::FT = FT(0.0)
    "Weight of specular2directional forward coefficient"
    sof::FT = FT(0.0)

    # dimension of nLevel
    "Probability of directly viewing a leaf in solar direction"
    Ps ::MArray{Tuple{nLevel},FT,1,nLevel} = MArray{Tuple{nLevel},FT}(undef)
    "Probability of directly viewing a leaf in viewing direction"
    Po ::MArray{Tuple{nLevel},FT,1,nLevel} = MArray{Tuple{nLevel},FT}(undef)
    "Bi-directional probability of directly viewing a leaf (solar->canopy->viewing)"
    Pso::MArray{Tuple{nLevel},FT,1,nLevel} = MArray{Tuple{nLevel},FT}(undef)

    # dimension of nIncl * nAzi
    "conversion factor fs to compute irradiance on inclined leaf"
    fs     ::MArray{Tuple{nIncl,nAzi},FT,2,nIncl_nAzi} = MArray{Tuple{nIncl,nAzi},FT}(undef)
    "abs(fs)"
    absfs  ::MArray{Tuple{nIncl,nAzi},FT,2,nIncl_nAzi} = MArray{Tuple{nIncl,nAzi},FT}(undef)
    "abs(fs*fo)"
    absfsfo::MArray{Tuple{nIncl,nAzi},FT,2,nIncl_nAzi} = MArray{Tuple{nIncl,nAzi},FT}(undef)
    "fs*fo"
    fsfo   ::MArray{Tuple{nIncl,nAzi},FT,2,nIncl_nAzi} = MArray{Tuple{nIncl,nAzi},FT}(undef)
    "conversion factor fo for angle towards observer (not sun like fs)"
    fo     ::MArray{Tuple{nIncl,nAzi},FT,2,nIncl_nAzi} = MArray{Tuple{nIncl,nAzi},FT}(undef)
    "Cosine of leaf azimuths"
    cosΘ_l ::MArray{Tuple{nIncl,nAzi},FT,2,nIncl_nAzi} = MArray{Tuple{nIncl,nAzi},FT}(undef)
    "cos of leaf azimuth sqared"
    cos2Θ_l::MArray{Tuple{nIncl,nAzi},FT,2,nIncl_nAzi} = MArray{Tuple{nIncl,nAzi},FT}(undef)

    # The following also depend on leaf reflectance and transmission.
    # Might go into a separate strcuture so that we can have it separately for thermal, SW and SIF?
    # dimension of nWL * nLayer
    "diffuse     backscatter scattering coefficient for diffuse  incidence"
    sigb::MArray{Tuple{nWL,nLayer},FT,2,nLayer_nWL} = MArray{Tuple{nWL,nLayer},FT}(undef)
    "diffuse     forward     scattering coefficient for diffuse  incidence"
    sigf::MArray{Tuple{nWL,nLayer},FT,2,nLayer_nWL} = MArray{Tuple{nWL,nLayer},FT}(undef)
    "diffuse     backscatter scattering coefficient for specular incidence"
    sb  ::MArray{Tuple{nWL,nLayer},FT,2,nLayer_nWL} = MArray{Tuple{nWL,nLayer},FT}(undef)
    "diffuse     forward     scattering coefficient for specular incidence"
    sf  ::MArray{Tuple{nWL,nLayer},FT,2,nLayer_nWL} = MArray{Tuple{nWL,nLayer},FT}(undef)
    "directional backscatter scattering coefficient for diffuse  incidence"
    vb  ::MArray{Tuple{nWL,nLayer},FT,2,nLayer_nWL} = MArray{Tuple{nWL,nLayer},FT}(undef)
    "directional forward     scattering coefficient for diffuse  incidence"
    vf  ::MArray{Tuple{nWL,nLayer},FT,2,nLayer_nWL} = MArray{Tuple{nWL,nLayer},FT}(undef)
    "bidirectional scattering coefficent (directional-directional)"
    w   ::MArray{Tuple{nWL,nLayer},FT,2,nLayer_nWL} = MArray{Tuple{nWL,nLayer},FT}(undef)
    "attenuation"
    a   ::MArray{Tuple{nWL,nLayer},FT,2,nLayer_nWL} = MArray{Tuple{nWL,nLayer},FT}(undef)
    "Effective layer transmittance (direct->diffuse)"
    Xsd ::MArray{Tuple{nWL,nLayer},FT,2,nLayer_nWL} = MArray{Tuple{nWL,nLayer},FT}(undef)
    "Effective layer transmittance (diffuse->diffuse)"
    Xdd ::MArray{Tuple{nWL,nLayer},FT,2,nLayer_nWL} = MArray{Tuple{nWL,nLayer},FT}(undef)

    # dimension of nWL * nLevel
    "Effective layer reflectance (direct->diffuse)"
    R_sd::MArray{Tuple{nWL,nLevel},FT,2,nLev_nWL} = MArray{Tuple{nWL,nLevel},FT}(undef)
    "Effective layer reflectance (diffuse->diffuse)"
    R_dd::MArray{Tuple{nWL,nLevel},FT,2,nLev_nWL} = MArray{Tuple{nWL,nLevel},FT}(undef)
    "Solar direct radiation per layer)"
    Es_ ::MArray{Tuple{nWL,nLevel},FT,2,nLev_nWL} = MArray{Tuple{nWL,nLevel},FT}(undef)
end




"""
    create_canopy_optical(FType, nWL::Int, nLayers::Int, nAzi::Int, nIncl::Int; using_marray=false)

Create a canopy optical properties struct with an `using_marray` option, given
- `FType` Floating number type
- `nWL` Number of wave length
- `nLayers` Number of canopy layers
- `nAzi` Number of arimuth angles
- `AIncl` Number of inclination angles
- `using_marray` If `true`, using MArray; else, using Array
"""
function create_canopy_optical(FType, nWL::Int, nLayers::Int, nAzi::Int, nIncl::Int; using_marray=false)
    if using_marray
        return CanopyOptiMArray{FType,nWL,nLayers,nLayers+1,nAzi,nIncl,nLayers*nWL,(nLayers+1)*nWL,nIncl*nAzi}()
    else
        return CanopyOptiArray{FType,nWL,nLayers,nLayers+1,nAzi,nIncl}()
    end
end








###############################################################################
#
# Canopy radiation
#
###############################################################################
"""
    struct CanopyRadiation{FT, nWl, nWlF, nIncl, nAzi, nLayers}

A struct for canopy radiation information

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct CanopyRadiation{FT, nWl, nWlF, nIncl, nAzi, nLayers}
    # Scalars
    "Integrated TOC outgoing flux `[W m⁻²]`"
    intEout            ::FT = FT(0.0)
    "Incident spectrally integrated total PAR `[mol m⁻² s⁻¹]`"
    incomingPAR        ::FT = FT(0.0)
    "Incident spectrally integrated direct PAR `[mol m⁻² s⁻¹]`"
    incomingPAR_direct ::FT = FT(0.0)
    "Incident spectrally integrated diffuse PAR `[mol m⁻² s⁻¹]`"
    incomingPAR_diffuse::FT = FT(0.0)
    "Net radiation of shaded soil `[W m⁻²]`"
    RnSoil_diffuse     ::FT = FT(0.0)
    "Net Short-wave radiation of sunlit soil `[W m⁻²]`"
    RnSoil_direct      ::FT = FT(0.0)
    "Net Short-wave radiation of soil (shaded + sunlit) `[W m⁻²]`"
    RnSoil             ::FT = FT(0.0)
    "Net long-wave radiation of soil (shaded + sunlit) `[W m⁻²]`"
    RnSoilLW           ::FT = FT(0.0)

    # Dim of nLayers
    "Net PAR of shaded leaves `[mol m⁻² s⁻¹]`"
    absPAR_shade   ::Array{FT,1} = zeros(FT, nLayers)
    "Net PAR by Cab+Car of shaded leaves `[moles m⁻² s⁻¹]`"
    absPAR_shadeCab::Array{FT,1} = zeros(FT, nLayers)
    "Spectrally integrated net absorbed direct radiation in each layer `[W m⁻²)]`"
    intNetSW_sunlit::Array{FT,1} = zeros(FT, nLayers)
    "Spectrally integrated net absorbed diffuse radiation in each layer `[W m⁻²)]`"
    intNetSW_shade ::Array{FT,1} = zeros(FT, nLayers)
    "Spectrally integrated net absorbed direct radiation in each layer `[W m⁻²)]`"
    intNetLW_sunlit::Array{FT,1} = zeros(FT, nLayers)
    "Spectrally integrated net absorbed diffuse radiation in each layer `[W m⁻²)]`"
    intNetLW_shade ::Array{FT,1} = zeros(FT, nLayers)
    "Leaf temperature (sunlit) `[K]`"
    T_sun          ::Array{FT,1} = zeros(FT, nLayers) .+ FT(298.15)
    "Leaf temperature (shaded) `[K]`"
    T_shade        ::Array{FT,1} = zeros(FT, nLayers) .+ FT(298.15)
    "Fluorescence yield for shaded leaves"
    ϕ_shade        ::Array{FT,1} =  ones(FT, nLayers) .* FT(0.01)
    "Sensible Heat flux H of shaded leaves `[W m⁻²]`"
    H_shade        ::Array{FT,1} = zeros(FT, nLayers)
    "Latent Heat flux LE of shaded leaves `[W m⁻²]`"
    LE_shade       ::Array{FT,1} = zeros(FT, nLayers)
    "NPQ of shaded leaves"
    NPQ_shade      ::Array{FT,1} = zeros(FT, nLayers)
    # TODO remove these?
    "GPP of shaded leaves `[μmol m⁻² s⁻¹]`"
    GPP_shade      ::Array{FT,1} = zeros(FT, nLayers)
    "gs of shaded leaves `[mol m⁻² s⁻¹]`"
    gs_shade       ::Array{FT,1} = zeros(FT, nLayers)
    "Leaf water potential of shaded leaves `[MPa]`"
    ψl_shade       ::Array{FT,1} = zeros(FT, nLayers)
    "Cc of shaded leaves `[µmol/mol]`"
    Cc_shade       ::Array{FT,1} = zeros(FT, nLayers)
    "internal CO₂ concentration of shaded leaves `[µmol/mol]`"
    Pi_shade       ::Array{FT,1} = zeros(FT, nLayers)

    # Dimension of wavelength
    "Short-wave TOC outgoing radiance in observation direction `[mW m⁻² nm⁻¹ sr⁻¹]`"
    Lo         ::Array{FT,1} = zeros(FT, nWl)
    "Short-wave TOC outgoing radiation `[mW m⁻² nm⁻¹]`"
    Eout       ::Array{FT,1} = zeros(FT, nWl)
    "Short-wave Albedo in viewing direction"
    alb_obs    ::Array{FT,1} = zeros(FT, nWl)
    "Short-wave Albedo for direct incoming radiation"
    alb_direct ::Array{FT,1} = zeros(FT, nWl)
    "Short-wave Albedo for diffuse incoming radiation"
    alb_diffuse::Array{FT,1} = zeros(FT, nWl)

    # Dimension of nLevel * nWavelengths
    "Upwelling diffuse short-wave radiation within canopy `[mW m⁻² nm⁻¹]`"
    E_up  ::Array{FT,2} = zeros(FT, (nWl,nLayers+1))
    "Downwelling diffuse short-wave radiation within canopy `[mW m⁻² nm⁻¹]`"
    E_down::Array{FT,2} = zeros(FT, (nWl,nLayers+1))

    # Dimension of nLayer * nWavelengths
    "Net absorbed direct radiation in each layer `[mW m⁻² nm⁻¹]`"
    netSW_sunlit   ::Array{FT,2} = zeros(FT, (nWl,nLayers))
    "net absorbed diffuse radiation in each layer `[mW m⁻² nm⁻¹]`"
    netSW_shade    ::Array{FT,2} = zeros(FT, (nWl,nLayers))


    # Dimension of nLeafInclination * nLeafAzimuth * nLayer
    "net PAR of sunlit leaves `[mol m⁻² s⁻¹]`"
    absPAR_sun   ::Array{FT,3} = zeros(FT, (nIncl,nAzi,nLayers))
    "net PAR by Cab+Car of sunlit leaves `[mol m⁻² s⁻¹]`"
    absPAR_sunCab::Array{FT,3} = zeros(FT, (nIncl,nAzi,nLayers))
    "Leaf temperature (sunlit) `[K]`"
    T_sun3D      ::Array{FT,3} = zeros(FT, (nIncl,nAzi,nLayers)) .+ FT(298.15)
    "Fluorescence yield for sunlit leaves"
    ϕ_sun        ::Array{FT,3} =  ones(FT, (nIncl,nAzi,nLayers)) .* FT(0.01)
    "Sensible Heat flux H of sunlit leaves `[W m⁻²]`"
    H_sun        ::Array{FT,3} = zeros(FT, (nIncl,nAzi,nLayers))
    "Latent Heat flux LE of sunlit leaves `[W m⁻²]`"
    LE_sun       ::Array{FT,3} = zeros(FT, (nIncl,nAzi,nLayers))
    "NPQ of sunlit leaves"
    NPQ_sun      ::Array{FT,3} = zeros(FT, (nIncl,nAzi,nLayers))
    # TODO remove these?
    "GPP of sunlit leaves `[μmol m⁻² s⁻¹]`"
    GPP_sun      ::Array{FT,3} = zeros(FT, (nIncl,nAzi,nLayers))
    "gs of sunlit leaves `[mol m⁻² s⁻¹]`"
    gs_sun       ::Array{FT,3} = zeros(FT, (nIncl,nAzi,nLayers))
    "Leaf water potential of sunlit leaves `[MPa]`"
    ψl_sun       ::Array{FT,3} = zeros(FT, (nIncl,nAzi,nLayers))
    "Cc of sunlit leaves `[µmol/mol]`"
    Cc_sun       ::Array{FT,3} = zeros(FT, (nIncl,nAzi,nLayers))
    "Internal CO₂ concentration of sunlit leaves `[µmol/mol]`"
    Pi_sun       ::Array{FT,3} = zeros(FT, (nIncl,nAzi,nLayers))

    # Fluorescence Output:
    "Hemispheric total outgoing SIF flux `[mW m⁻² nm⁻¹]`)"
    SIF_hemi         ::Array{FT,1} = zeros(FT, nWlF)
    "Observer-direction outgoing SIF radiance  (mW m⁻² nm⁻¹ sr⁻¹))"
    SIF_obs          ::Array{FT,1} = zeros(FT, nWlF)
    "Observer-direction outgoing SIF radiance, sunlit leaves  (mW m⁻² nm⁻¹ sr⁻¹))"
    SIF_obs_sunlit   ::Array{FT,1} = zeros(FT, nWlF)
    "Observer-direction outgoing SIF radiance, shaded leaves  (mW m⁻² nm⁻¹ sr⁻¹))"
    SIF_obs_shaded   ::Array{FT,1} = zeros(FT, nWlF)
    "Observer-direction outgoing SIF radiance, scattered   (mW m⁻² nm⁻¹ sr⁻¹))"
    SIF_obs_scattered::Array{FT,1} = zeros(FT, nWlF)
    "Observer-direction outgoing SIF radiance, soil-reflected  (mW m⁻² nm⁻¹ sr⁻¹))"
    SIF_obs_soil     ::Array{FT,1} = zeros(FT, nWlF)
    "Total SIF sum of layer sources  `[mW m⁻² nm⁻¹]`)"
    SIF_sum          ::Array{FT,1} = zeros(FT, nWlF)
end








###############################################################################
#
# Incoming radiation
#
###############################################################################
"""
    AbstractIncomingRadiation

Abstract leaf biological parameters type
"""
abstract type AbstractIncomingRadiation end




"""
    struct IncomingRadiationArray{FT}

Incoming radiation information.

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct IncomingRadiationArray{FT} <: AbstractIncomingRadiation
    "Wavelength `[nm]`"
    wl       ::Array{FT,1}
    "Direct incoming radiation `[mW m⁻² nm⁻¹]`"
    E_direct ::Array{FT,1}
    "Diffuse incoming radiation `[mW m⁻² nm⁻¹]`"
    E_diffuse::Array{FT,1}
end




"""
    struct IncomingRadiationMArray{FT, N}

Incoming radiation information.

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct IncomingRadiationMArray{FT, N} <: AbstractIncomingRadiation
    "Wavelength `[nm]`"
    wl       ::MArray{Tuple{N},FT,1,N}
    "Direct incoming radiation `[mW m⁻² nm⁻¹]`"
    E_direct ::MArray{Tuple{N},FT,1,N}
    "Diffuse incoming radiation `[mW m⁻² nm⁻¹]`"
    E_diffuse::MArray{Tuple{N},FT,1,N}
end




"""
    create_incoming_radiation(FType, swl:Array, file; using_marray=false)

Create an `AbstractIncomingRadiation` struct, given
- `FType` Floating number type
- `swl` Standard wave length
- `file` Input file name
- `using_marray` If `true`, use MArray; else, using Array
"""
function create_incoming_radiation(FType, swl::Array, file::String = file_Sun; using_marray::Bool=false)
    N = length(swl)-1

    # Read data
    _suni  = matread(file)["sun"]
    _wl    = _suni["wl"      ]
    _Edir  = _suni["Edirect" ]
    _Ediff = _suni["Ediffuse"]

    # create arrays
    if using_marray
        wl    = MArray{Tuple{N},FType}(undef)
        Edir  = MArray{Tuple{N},FType}(undef)
        Ediff = MArray{Tuple{N},FType}(undef)
    else
        wl    = zeros(FType, N)
        Edir  = zeros(FType, N)
        Ediff = zeros(FType, N)
    end

    # fill in the arrays
    # println("Reading Optical Parameters from ", swl[1], " to ", swl[end], " length: ", length(swl))
    for i in 1:N
        wo = findall( (_wl.>=swl[i]) .& (_wl.<swl[i+1]) )
        if length(wo)==0
            println("Warning, some wavelengths out of bounds ", swl[i])
        end
        wl[i]    = mean(   _wl[wo])
        Edir[i]  = mean( _Edir[wo])
        Ediff[i] = mean(_Ediff[wo])
    end

    # create struct from the arrays
    if using_marray
        return IncomingRadiationMArray{FType,N}(wl, Edir, Ediff)
    else
        return IncomingRadiationArray{FType}(wl, Edir, Ediff)
    end
end








###############################################################################
#
# Leaf biological parameters
#
###############################################################################
"""
    AbstractLeafBio

Abstract leaf biological parameters type
"""
abstract type AbstractLeafBio end




"""
    struct LeafBioArray{FT, nWl, nWle, nWlf}

A struct of leaf biological parameters using Array

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct LeafBioArray{FT, nWl, nWle, nWlf} <: AbstractLeafBio
    "Leaf structure parameter"
    N   ::FT = FT(1.4  )
    "Chlorophyll a+b content `[µg cm⁻²]`"
    Cab ::FT = FT(40.0 )
    "Carotenoid content `[µg cm⁻²]`"
    Car ::FT = FT(10.0 )
    "Anthocynanin content `[µg cm⁻²]`"
    Ant ::FT = FT(0.0  )
    "Senescent material fraction"
    Cs  ::FT = FT(0.0  )
    "Equivalent water thickness `[cm]`"
    Cw  ::FT = FT(0.009)
    "Dry matter content (dry leaf mass per unit area) `[g cm⁻²]`"
    Cm  ::FT = FT(0.012)
    "Fractionation between Zeaxanthin and Violaxanthin in Car (1=all Zeaxanthin) (-)"
    Cx  ::FT = FT(0.0  )
    "Broadband thermal reflectance (-)"
    ρ_LW::FT = FT(0.01 )
    "Broadband thermal transmission (-)"
    τ_LW::FT = FT(0.01 )
    "Leaf fluorescence efficiency (Fo standard)"
    fqe ::FT = FT(1.0  )

    "Shortwave leaf reflectance"
    ρ_SW       ::Array{FT,1} = zeros(FT, nWl)
    "Shortwave leaf transmission"
    τ_SW       ::Array{FT,1} = zeros(FT, nWl)
    "Relative absorbtion by Chlorophyll+Car"
    kChlrel    ::Array{FT,1} = zeros(FT, nWl)
    "Relative absorbtion by Chlorophyll"
    kChlrel_old::Array{FT,1} = zeros(FT, nWl)
    "Fluorescence excitation matrix backwards"
    Mb         ::Array{FT,2} = zeros(FT,(nWlf,nWle))
    "Fluorescence excitation matrix forwards"
    Mf         ::Array{FT,2} = zeros(FT,(nWlf,nWle))
    "Doubling adding layers"
    ndub      ::Int = 10
end




"""
    struct LeafBioMArray{FT, nWl, nWle, nWlf, nWLe_nWLf}

A struct of leaf biological parameters using Array

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct LeafBioMArray{FT, nWl, nWle, nWlf, nWLe_nWLf} <: AbstractLeafBio
    "Leaf structure parameter"
    N   ::FT = FT(1.4  )
    "Chlorophyll a+b content `[µg cm⁻²]`"
    Cab ::FT = FT(40.0 )
    "Carotenoid content `[µg cm⁻²]`"
    Car ::FT = FT(10.0 )
    "Anthocynanin content `[µg cm⁻²]`"
    Ant ::FT = FT(0.0  )
    "Senescent material fraction"
    Cs  ::FT = FT(0.0  )
    "Equivalent water thickness `[cm]`"
    Cw  ::FT = FT(0.009)
    "Dry matter content (dry leaf mass per unit area) `[g cm⁻²]`"
    Cm  ::FT = FT(0.012)
    "Fractionation between Zeaxanthin and Violaxanthin in Car (1=all Zeaxanthin)"
    Cx  ::FT = FT(0.0  )
    "Broadband thermal reflectance"
    ρ_LW::FT = FT(0.01 )
    "Broadband thermal transmission"
    τ_LW::FT = FT(0.01 )
    "Leaf fluorescence efficiency (Fo standard)"
    fqe ::FT = FT(1.0  )

    "Shortwave leaf reflectance"
    ρ_SW       ::MArray{Tuple{nWl},FT,1,nWl} = MArray{Tuple{nWl},FT}(undef)
    "Shortwave leaf transmission"
    τ_SW       ::MArray{Tuple{nWl},FT,1,nWl} = MArray{Tuple{nWl},FT}(undef)
    "Relative absorbtion by Chlorophyll+Car"
    kChlrel    ::MArray{Tuple{nWl},FT,1,nWl} = MArray{Tuple{nWl},FT}(undef)
    "Relative absorbtion by Chlorophyll"
    kChlrel_old::MArray{Tuple{nWl},FT,1,nWl} = MArray{Tuple{nWl},FT}(undef)
    "Fluorescence backward excitation matrix "
    Mb::MArray{Tuple{nWlf,nWle},FT,2,nWLe_nWLf} = MArray{Tuple{nWlf,nWle},FT}(undef)
    "Fluorescence forward excitation matrix "
    Mf::MArray{Tuple{nWlf,nWle},FT,2,nWLe_nWLf} = MArray{Tuple{nWlf,nWle},FT}(undef)
    "Doubling adding layers"
    ndub      ::Int = 10
end




"""
    create_leaf_bio(FT, nWl::Int, nWle::Int, nWlf::Int; using_marray=false)

Create a leaf biological parameters struct with an `using_marray` option, given
- `FType` Floating number type
- `nWl` Number of wave length
- `nWle` Number of excitation wave length
- `nWlf` Number of fluorescence wave length
- `using_marray` If `true`, using MArray; else, using Array

Returns a [`LeafBioMArray`](@ref) or [`LeafBioArray`](@ref) type struct.
"""
function create_leaf_bio(FType, nWl::Int, nWle::Int, nWlf::Int; using_marray=false)
    if using_marray
        return LeafBioMArray{FType, nWl, nWle, nWlf, nWle*nWlf}()
    else
        return LeafBioArray{FType, nWl, nWle, nWlf}()
    end
end








###############################################################################
#
# Leaf optical parameters
#
###############################################################################
"""
AbstractLeafOptiPara

Abstract leaf optical parameter
"""
abstract type AbstractLeafOptiPara end




"""
    struct LeafOptiParaArray{FT, N}

Struct for leaf optical properties using Array

# Fields
$(DocStringExtensions.FIELDS)
"""
struct LeafOptiParaArray{FT} <: AbstractLeafOptiPara
    # TODO Add explanations to each field
    nr    ::Array{FT,1}
    Km    ::Array{FT,1}
    Kab   ::Array{FT,1}
    Kant  ::Array{FT,1}
    Kcar  ::Array{FT,1}
    Kw    ::Array{FT,1}
    KBrown::Array{FT,1}
    phi   ::Array{FT,1}
    KcaV  ::Array{FT,1}
    KcaZ  ::Array{FT,1}
    lambda::Array{FT,1}
end




"""
    struct LeafOptiParaMArray{FT, N}

Struct for leaf optical properties using MArray

# Fields
$(DocStringExtensions.FIELDS)
"""
struct LeafOptiParaMArray{FT, N} <: AbstractLeafOptiPara
    # TODO Add explanations to each field
    nr    ::MArray{Tuple{N},FT,1,N}
    Km    ::MArray{Tuple{N},FT,1,N}
    Kab   ::MArray{Tuple{N},FT,1,N}
    Kant  ::MArray{Tuple{N},FT,1,N}
    Kcar  ::MArray{Tuple{N},FT,1,N}
    Kw    ::MArray{Tuple{N},FT,1,N}
    KBrown::MArray{Tuple{N},FT,1,N}
    phi   ::MArray{Tuple{N},FT,1,N}
    KcaV  ::MArray{Tuple{N},FT,1,N}
    KcaZ  ::MArray{Tuple{N},FT,1,N}
    lambda::MArray{Tuple{N},FT,1,N}
end




"""
    create_opti_par(FType, swl::Array, file; using_marray=false)

Create an `AbstractLeafOptiPara` struct, given
- `FType` Floating number type
- `swl` Standard wave length
- `file` Input file name
- `using_marray` If `true`, use MArray; else, using Array
"""
function create_opti_par(FType, swl::Array, file::String=file_Opti; using_marray=false)
    N = length(swl)-1

    # reading data
    _opti   = matread(file)["optipar"]
    _nr     = _opti["nr"  ]
    _Km     = _opti["Kdm" ]
    _Kab    = _opti["Kab" ]
    _Kant   = _opti["Kant"]
    _Kcar   = _opti["Kca" ]
    _Kw     = _opti["Kw"  ]
    _KBrown = _opti["Ks"  ]
    _phi    = _opti["phi" ]
    _KcaV   = _opti["KcaV"]
    _KcaZ   = _opti["KcaZ"]
    _lambda = _opti["wl"  ]

    # create data to parse
    if using_marray
        nr     = MArray{Tuple{N},FType}(undef)
        Km     = MArray{Tuple{N},FType}(undef)
        Kab    = MArray{Tuple{N},FType}(undef)
        Kant   = MArray{Tuple{N},FType}(undef)
        Kcar   = MArray{Tuple{N},FType}(undef)
        Kw     = MArray{Tuple{N},FType}(undef)
        KBrown = MArray{Tuple{N},FType}(undef)
        phi    = MArray{Tuple{N},FType}(undef)
        KcaV   = MArray{Tuple{N},FType}(undef)
        KcaZ   = MArray{Tuple{N},FType}(undef)
        lambda = MArray{Tuple{N},FType}(undef)
    else
        nr     = zeros(FType, N)
        Km     = zeros(FType, N)
        Kab    = zeros(FType, N)
        Kant   = zeros(FType, N)
        Kcar   = zeros(FType, N)
        Kw     = zeros(FType, N)
        KBrown = zeros(FType, N)
        phi    = zeros(FType, N)
        KcaV   = zeros(FType, N)
        KcaZ   = zeros(FType, N)
        lambda = zeros(FType, N)
    end

    # fill in the data arrays
    # println("Reading Optical Parameters from ", swl[1], " to ", swl[end], " length: ", length(swl))
    @inbounds for i in 1:N
        wo = findall( (_lambda.>=swl[i]) .& (_lambda.<swl[i+1]) )
        if length(wo)==0
            println("Warning, some wavelengths out of bounds ", swl[i])
        end

        nr[i]     = mean(    _nr[wo])
        Km[i]     = mean(    _Km[wo])
        Kab[i]    = mean(   _Kab[wo])
        Kant[i]   = mean(  _Kant[wo])
        Kcar[i]   = mean(  _Kcar[wo])
        Kw[i]     = mean(    _Kw[wo])
        KBrown[i] = mean(_KBrown[wo])
        phi[i]    = mean(   _phi[wo])
        KcaV[i]   = mean(  _KcaV[wo])
        KcaZ[i]   = mean(  _KcaZ[wo])
        lambda[i] = mean(_lambda[wo])
    end

    # return the created struct
    if using_marray
        return LeafOptiParaMArray{FType,N}(nr, Km, Kab, Kant, Kcar, Kw, KBrown, phi, KcaV, KcaZ, lambda)
    else
        return LeafOptiParaArray{FType}(nr, Km, Kab, Kant, Kcar, Kw, KBrown, phi, KcaV, KcaZ, lambda)
    end
end








###############################################################################
#
# Soil optical parameters
#
###############################################################################
"""
    struct SoilOpti{FT}

A struct of soil optical parameters

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct SoilOpti{FT}
    "Wavelength `[nm]`"
    wl        ::Array{FT,1}
    "Shortwave albedo"
    albedo_SW ::Array{FT,1}
    "Longwave albedo"
    albedo_LW ::Array{FT,1}
    "Soil surface temperature `[K]`"
    soil_skinT::FT
end








###############################################################################
#
# Solar angle type
#
###############################################################################
"""
    struct SolarAngles{FT}

Struct for observation and solar angles

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct SolarAngles{FT}
    "Solar Zenith Angle `[degree]`"
    tts::FT = FT(30.0)
    "Viewing Zenith Angle in `[degree]`"
    tto::FT = FT(0.0 )
    "relative azimuth in `[degree]`"
    psi::FT = FT(0.0 )
end








###############################################################################
#
# Wave length parameter set
#
###############################################################################
"""
    AbstractWLParaSet

An abstract type for wave length parameter set, which corresponds to the AbstractLeafOptiPara
"""
abstract type AbstractWLParaSet end




"""
    struct WLParaSetArray{FT}

Struct for pre-set wave length parameters.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct WLParaSetArray{FT} <: AbstractWLParaSet
    # Wave length (WL) boundaries
    "Minimal WL for PAR `[nm]`"
    minwlPAR::FT = FT(400.0)
    "Maximal WL for PAR `[nm]`"
    maxwlPAR::FT = FT(700.0)
    "Minimal WL for SIF excitation `[nm]`"
    minwle  ::FT = FT(400.0)
    "Maximal WL for SIF excitation `[nm]`"
    maxwle  ::FT = FT(750.0)
    "Minimal WL for SIF emission/fluorescence `[nm]`"
    minwlf  ::FT = FT(640.0)
    "Maximal WL for SIF emission/fluorescence `[nm]` "
    maxwlf  ::FT = FT(850.0)

    # Wave length lists
    "Standard wave length `[nm]`"
    swl::Array{FT,1} = [collect(FT(400.0):FT(10.0):FT( 650.1));
                        collect(FT(655.0):FT( 5.0):FT( 770.1));
                        collect(FT(780.0):FT(25.0):FT(2400.1))]
    "Differential wavelength"
    dwl::Array{FT,1} = diff(swl)

    "Leaf optical parameter set"
    optis::LeafOptiParaArray = create_opti_par(FT, swl, file_Opti; using_marray=false)

    "Wave length `[nm]`"
    wl::Array{FT,1} = optis.lambda

    "Length of wl"
    nwl ::Int         = length(wl)
    "Index of wle in wl"
    Iwle::Array       = findall( (wl .>= minwle) .& (wl .<= maxwle) )
    "Index of wlf in wl"
    Iwlf::Array       = findall( (wl .>= minwlf) .& (wl .<= maxwlf) )
    "length of wle"
    nWlE::Int         = length(Iwle)
    "length of wlf"
    nWlF::Int         = length(Iwlf)
    "index of wlPAR in wl"
    iPAR::Array       = findall( (wl .>= minwlPAR) .& (wl .<= maxwlPAR) )
    "excitation wave length `[nm]`"
    wle ::Array{FT,1} = wl[Iwle]
    "Fluorescence wave length `[nm]`"
    wlf ::Array{FT,1} = wl[Iwlf]
end




"""
    struct WLParaSetMArray{FT}

Struct for pre-set wave length parameters.

# TODO add MArray version?

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct WLParaSetMArray{FT} <: AbstractWLParaSet
    # Wave length (WL) boundaries
    "Minimal WL for PAR `[nm]`"
    minwlPAR::FT = FT(400.0)
    "Maximal WL for PAR `[nm]`"
    maxwlPAR::FT = FT(700.0)
    "Minimal WL for SIF excitation `[nm]`"
    minwle  ::FT = FT(400.0)
    "Maximal WL for SIF excitation `[nm]`"
    maxwle  ::FT = FT(750.0)
    "Minimal WL for SIF emission/fluorescence `[nm]`"
    minwlf  ::FT = FT(640.0)
    "Maximal WL for SIF emission/fluorescence `[nm]` "
    maxwlf  ::FT = FT(850.0)

    # Wave length lists
    "Standard wave length `[nm]`"
    swl::Array{FT,1} = [collect(FT(400.0):FT(10.0):FT( 650.1));
                        collect(FT(655.0):FT( 5.0):FT( 770.1));
                        collect(FT(780.0):FT(25.0):FT(2400.1))]
    "Differential wavelength"
    dwl::Array{FT,1} = diff(swl)

    "Leaf optical parameter set"
    optis::LeafOptiParaMArray = create_opti_par(FT, swl, file_Opti; using_marray=true)

    "Wave length `[nm]`"
    wl::MArray = optis.lambda

    "Length of wl"
    nwl ::Int         = length(wl)
    "Index of wle in wl"
    Iwle::Array       = findall( (wl .>= minwle) .& (wl .<= maxwle) )
    "Index of wlf in wl"
    Iwlf::Array       = findall( (wl .>= minwlf) .& (wl .<= maxwlf) )
    "length of wle"
    nWlE::Int         = length(Iwle)
    "length of wlf"
    nWlF::Int         = length(Iwlf)
    "index of wlPAR in wl"
    iPAR::Array       = findall( (wl .>= minwlPAR) .& (wl .<= maxwlPAR) )
    "excitation wave length `[nm]`"
    wle ::Array{FT,1} = wl[Iwle]
    "Fluorescence wave length `[nm]`"
    wlf ::Array{FT,1} = wl[Iwlf]
end





"""
    create_wl_para_set(FType; using_marray::Bool=false)

Create a pre-set struct of wave length settings, given
- `Ftype` Floating number type
- `using_marray` If true, use MArray; else, use Array.

# TODO add MArray version?
"""
function create_wl_para_set(FType; using_marray::Bool=false)
    if using_marray
        return WLParaSetMArray{FType}()
    else
        return WLParaSetArray{FType}()
    end
end
