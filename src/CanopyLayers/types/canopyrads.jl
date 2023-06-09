###############################################################################
#
# Canopy radiation
#
###############################################################################
"""
    mutable struct CanopyRads{FT}

A struct for canopy radiation information

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct CanopyRads{FT}
    # local storage of dimension information
    "Number of azimuth angles"
    nAzi::Int = 36
    "Number of inclination agles"
    nIncl::Int = 9
    "Number of canopy layers"
    nLayer::Int = 5
    "Number of canopy levels"
    nLevel::Int = nLayer+1
    "Number of wave lengths"
    nWL::Int = 10
    "Number of wave lengths for SIF"
    nWLF::Int = 10

    # Scalars
    "Integrated TOC outgoing flux `[W m⁻²]`"
    intEout::FT = 0
    "Incident spectrally integrated total PAR `[mol m⁻² s⁻¹]`"
    incomingPAR::FT = 0
    "Incident spectrally integrated direct PAR `[mol m⁻² s⁻¹]`"
    incomingPAR_direct::FT = 0
    "Incident spectrally integrated diffuse PAR `[mol m⁻² s⁻¹]`"
    incomingPAR_diffuse::FT = 0
    "Net radiation of shaded soil `[W m⁻²]`"
    RnSoil_diffuse::FT = 0
    "Net Short-wave radiation of sunlit soil `[W m⁻²]`"
    RnSoil_direct::FT = 0
    "Net Short-wave radiation of soil (shaded + sunlit) `[W m⁻²]`"
    RnSoil::FT = 0
    "Net long-wave radiation of soil (shaded + sunlit) `[W m⁻²]`"
    RnSoilLW::FT = 0

    # Dim of nLayer
    "Net PAR of shaded leaves `[mol m⁻² s⁻¹]`"
    absPAR_shade::Vector{FT} = zeros(FT, nLayer)
    "Net PAR by Cab+Car of shaded leaves `[moles m⁻² s⁻¹]`"
    absPAR_shadeCab::Vector{FT} = zeros(FT, nLayer)
    "Spectrally integrated net absorbed direct radiation in each layer per leaf area `[W m⁻²]`"
    intNetSW_sunlit::Vector{FT} = zeros(FT, nLayer)
    "Spectrally integrated net absorbed diffuse radiation in each layer per leaf area `[W m⁻²]`"
    intNetSW_shade::Vector{FT} = zeros(FT, nLayer)
    "Spectrally integrated net absorbed direct radiation in each layer per leaf area `[W m⁻²]`"
    intNetLW_sunlit::Vector{FT} = zeros(FT, nLayer)
    "Spectrally integrated net absorbed diffuse radiation in each layer per leaf area `[W m⁻²]`"
    intNetLW_shade::Vector{FT} = zeros(FT, nLayer)
    "Leaf temperature (sunlit) `[K]`"
    T_sun::Vector{FT} = zeros(FT, nLayer) .+ FT(298.15)
    "Leaf temperature (shaded) `[K]`"
    T_shade::Vector{FT} = zeros(FT, nLayer) .+ FT(298.15)
    "Fluorescence yield for shaded leaves"
    ϕ_shade::Vector{FT} =  ones(FT, nLayer) .* FT(0.01)
    "Sensible Heat flux H of shaded leaves `[W m⁻²]`"
    H_shade::Vector{FT} = zeros(FT, nLayer)
    "Latent Heat flux LE of shaded leaves `[W m⁻²]`"
    LE_shade::Vector{FT} = zeros(FT, nLayer)
    "NPQ of shaded leaves"
    NPQ_shade::Vector{FT} = zeros(FT, nLayer)
    # TODO remove these?
    "GPP of shaded leaves `[μmol m⁻² s⁻¹]`"
    GPP_shade::Vector{FT} = zeros(FT, nLayer)
    "gs of shaded leaves `[mol m⁻² s⁻¹]`"
    gs_shade::Vector{FT} = zeros(FT, nLayer)
    "Leaf water potential of shaded leaves `[MPa]`"
    ψl_shade::Vector{FT} = zeros(FT, nLayer)
    "Cc of shaded leaves `[µmol/mol]`"
    Cc_shade::Vector{FT} = zeros(FT, nLayer)
    "internal CO₂ concentration of shaded leaves `[µmol/mol]`"
    Pi_shade::Vector{FT} = zeros(FT, nLayer)

    # Dimension of wavelength
    "Short-wave TOC outgoing radiance in observation direction `[mW m⁻² nm⁻¹ sr⁻¹]`"
    Lo::Vector{FT} = zeros(FT, nWL)
    "Short-wave TOC outgoing radiation `[mW m⁻² nm⁻¹]`"
    Eout::Vector{FT} = zeros(FT, nWL)
    "Short-wave Albedo in viewing direction"
    alb_obs::Vector{FT} = zeros(FT, nWL)
    "Short-wave Albedo for direct incoming radiation"
    alb_direct::Vector{FT} = zeros(FT, nWL)
    "Short-wave Albedo for diffuse incoming radiation"
    alb_diffuse::Vector{FT} = zeros(FT, nWL)

    # Dimension of nLevel * nWavelengths
    "Upwelling diffuse short-wave radiation within canopy `[mW m⁻² nm⁻¹]`"
    E_up::Matrix{FT} = zeros(FT, (nWL,nLevel))
    "Downwelling diffuse short-wave radiation within canopy `[mW m⁻² nm⁻¹]`"
    E_down::Matrix{FT} = zeros(FT, (nWL,nLevel))

    # Dimension of nLayer * nWavelengths
    "Net absorbed direct radiation in each layer `[mW m⁻² nm⁻¹]`"
    netSW_sunlit::Matrix{FT} = zeros(FT, (nWL,nLayer))
    "net absorbed diffuse radiation in each layer `[mW m⁻² nm⁻¹]`"
    netSW_shade::Matrix{FT} = zeros(FT, (nWL,nLayer))

    # Dimension of nLeafInclination * nLeafAzimuth * nLayer
    "net PAR of sunlit leaves `[mol m⁻² s⁻¹]`"
    absPAR_sun::Array{FT,3} = zeros(FT, (nIncl,nAzi,nLayer))
    "net PAR by Cab+Car of sunlit leaves `[mol m⁻² s⁻¹]`"
    absPAR_sunCab::Array{FT,3} = zeros(FT, (nIncl,nAzi,nLayer))
    "Leaf temperature (sunlit) `[K]`"
    T_sun3D::Array{FT,3} = zeros(FT, (nIncl,nAzi,nLayer)) .+ FT(298.15)
    "Fluorescence yield for sunlit leaves"
    ϕ_sun::Array{FT,3} =  ones(FT, (nIncl,nAzi,nLayer)) .* FT(0.01)
    "Sensible Heat flux H of sunlit leaves `[W m⁻²]`"
    H_sun::Array{FT,3} = zeros(FT, (nIncl,nAzi,nLayer))
    "Latent Heat flux LE of sunlit leaves `[W m⁻²]`"
    LE_sun::Array{FT,3} = zeros(FT, (nIncl,nAzi,nLayer))
    "NPQ of sunlit leaves"
    NPQ_sun::Array{FT,3} = zeros(FT, (nIncl,nAzi,nLayer))
    # TODO remove these?
    "GPP of sunlit leaves `[μmol m⁻² s⁻¹]`"
    GPP_sun::Array{FT,3} = zeros(FT, (nIncl,nAzi,nLayer))
    "gs of sunlit leaves `[mol m⁻² s⁻¹]`"
    gs_sun::Array{FT,3} = zeros(FT, (nIncl,nAzi,nLayer))
    "Leaf water potential of sunlit leaves `[MPa]`"
    ψl_sun::Array{FT,3} = zeros(FT, (nIncl,nAzi,nLayer))
    "Cc of sunlit leaves `[µmol/mol]`"
    Cc_sun::Array{FT,3} = zeros(FT, (nIncl,nAzi,nLayer))
    "Internal CO₂ concentration of sunlit leaves `[µmol/mol]`"
    Pi_sun::Array{FT,3} = zeros(FT, (nIncl,nAzi,nLayer))

    # Fluorescence Output:
    "Hemispheric total outgoing SIF flux `[mW m⁻² nm⁻¹]`)"
    SIF_hemi::Vector{FT} = zeros(FT, nWLF)
    "Observer-direction outgoing SIF radiance (mW m⁻² nm⁻¹ sr⁻¹))"
    SIF_obs::Vector{FT} = zeros(FT, nWLF)
    "Observer-direction outgoing SIF radiance, sunlit leaves (mW m⁻² nm⁻¹ sr⁻¹)"
    SIF_obs_sunlit::Vector{FT} = zeros(FT, nWLF)
    "Observer-direction outgoing SIF radiance, shaded leaves (mW m⁻² nm⁻¹ sr⁻¹)"
    SIF_obs_shaded::Vector{FT} = zeros(FT, nWLF)
    "Observer-direction outgoing SIF radiance, scattered (mW m⁻² nm⁻¹ sr⁻¹)"
    SIF_obs_scattered::Vector{FT} = zeros(FT, nWLF)
    "Observer-direction outgoing SIF radiance, soil-reflected (mW m⁻² nm⁻¹ sr⁻¹)"
    SIF_obs_soil::Vector{FT} = zeros(FT, nWLF)
    "Total SIF sum of layer sources  `[mW m⁻² nm⁻¹]`)"
    SIF_sum::Vector{FT} = zeros(FT, nWLF)
end
