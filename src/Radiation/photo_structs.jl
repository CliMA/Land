using Parameters
using DocStringExtensions
using StaticArrays

# Fixed path right now here
file_Opti = joinpath(@__DIR__, "Optipar2017_ProspectD.mat")
file_Sun = joinpath(@__DIR__, "sun.mat")

# Struct for observation and solar angles
"""
    Struct for observation and solar angles

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct struct_angles{FT<:Number}
    "Solar Zenith Angle in degrees"
    tts::FT=30
    "Viewing Zenith Angle in degrees"
    tto::FT=0
    "relative azimuth in degrees"
    psi::FT=0
end

# Struct for leaf optical properties
struct optipar{FT<:AbstractFloat,N}
    nr::MArray{Tuple{N}, FT,1,N}
    Km::MArray{Tuple{N}, FT,1,N}
    Kab::MArray{Tuple{N}, FT,1,N}
    Kant::MArray{Tuple{N}, FT,1,N}
    Kcar::MArray{Tuple{N}, FT,1,N}
    Kw::MArray{Tuple{N}, FT,1,N}
    KBrown::MArray{Tuple{N}, FT,1,N}
    phi::MArray{Tuple{N}, FT,1,N}
    KcaV::MArray{Tuple{N}, FT,1,N}
    KcaZ::MArray{Tuple{N}, FT,1,N}
    lambda::MArray{Tuple{N}, FT,1,N}
end

"""
    incomingRadiation

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct incomingRadiation{FT<:AbstractFloat}
    "Wavelength (nm)"
    wl::Array{FT,1}
    "Direct incoming radiation `[mW m^-2 μm^-1]`"
    E_direct::Array{FT,1}
    "Diffuse incoming radiation `[mW m^-2 μm^-1]`"
    E_diffuse::Array{FT,1}

end

"""
    leafbio

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct leafbio{FT<:AbstractFloat,nWl,nWle, nWlf,nWLe_nWLf}
    "Leaf structure parameter (-)"
    N::FT    = 1.4      
    "Chlorophyll a+b content (µg cm-2)"
    Cab::FT  = 40.0      
    "Carotenoid content (µg cm-2)"
    Car::FT  = 10.0      
    "Anthocynanin content (µg cm-2)"
    Ant::FT  = 0.0       
    "Senescent material fraction (-)"
    Cs::FT   = 0.0      
    "Equivalent water thickness (cm)"
    Cw::FT   = 0.009     
    "Dry matter content (dry leaf mass per unit area) (g cm-2)"
    Cm::FT   = 0.012      
    "Fractionation between Zeaxanthin and Violaxanthin in Car (1=all Zeaxanthin) (-)"
    Cx::FT   = 0.0       
    "Broadband thermal reflectance (-)"
    ρ_LW::FT = 0.01      
    "Broadband thermal transmission (-)"
    τ_LW::FT = 0.01      
    "Leaf fluorescence efficiency (Fo standard)"
    fqe::FT = 1       

    "shortwave leaf reflectance (-)"
    ρ_SW::MArray{Tuple{nWl}, FT,1,nWl}           = MArray{Tuple{nWl}, FT}(undef);    # | -          | (0.0, 1.0)  | "shortwave reflectance"
    "shortwave leaf transmission (-)"
    τ_SW::MArray{Tuple{nWl}, FT,1,nWl}           = MArray{Tuple{nWl}, FT}(undef);   # | -          | (0.0, 1.0)  | "shortwave transmission"
    "relative absorbtion by Chlorophyll+Car (-)"
    kChlrel::MArray{Tuple{nWl}, FT,1,nWl}        = MArray{Tuple{nWl}, FT}(undef); # | -          | (0.0, 1.0)  | "relative absorbtion by Chlorophyll"
    "relative absorbtion by Chlorophyll (-)"
    kChlrel_old::MArray{Tuple{nWl}, FT,1,nWl}    = MArray{Tuple{nWl}, FT}(undef); # | -          | (0.0, 1.0)  | "relative absorbtion by Chlorophyll"
    "Fluorescence excitation matrix backwards (-)"
    Mb::MArray{Tuple{nWlf,nWle}, FT,2,nWLe_nWLf} = MArray{Tuple{nWlf,nWle}, FT}(undef);     # | -          | (0.0, 1.0)  | "Fluorescence excitation matrix backwards"
    "Fluorescence excitation matrix forwards (-)"
    Mf::MArray{Tuple{nWlf,nWle}, FT,2,nWLe_nWLf} = MArray{Tuple{nWlf,nWle}, FT}(undef);     # | -          | (0.0, 1.0)  | "Fluorescence excitation matrix forwards"
end

"""
    struct_soil

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct struct_soil{FT<:AbstractFloat}
    "Wavelength (nm)"
    wl::Array{FT,1}
    "shortwave albedo"
    albedo_SW::Array{FT,1}    # | -          | (0.0, 1.0)  | "shortwave albedo"
    "longwave albedo"
    albedo_LW::Array{FT,1}    # | -          | (0.0, 1.0)  | "longwave albedo"
    "Soil surface temperature (K)"
    soil_skinT::FT    # | -          | (0.0, 1.0)  | "longwave albedo"
end

"""
    struct_canopyRadiation

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct struct_canopyRadiation{FT<:AbstractFloat,nWl,nWlF,nIncl,nAzi,nLayers}
    # Scalars
    "integrated TOC outgoing flux `[W m^-2]`"
    intEout::FT = 0;                       # | W m^-2               | (0.0, 2500)  | "integrated TOC outgoing flux"
    "incident spectrally integrated total PAR `[moles m^-2 s^-1]`"
    incomingPAR::FT = 0;                           # | moles m^-2 s^-1      | (0.0, 2.5e-9)| "incident spectrally integrated total PAR"
    "incident spectrally integrated direct PAR `[moles m^-2 s^-1]`"
    incomingPAR_direct::FT = 0;                    # | moles m^-2 s^-1      | (0.0, 2.5e-9)| "incident spectrally integrated direct PAR"
    "incident spectrally integrated diffuse PAR `[moles m^-2 s^-1]`"
    incomingPAR_diffuse::FT = 0;                  # | moles m^-2 s^-1      | (0.0, 2.5e-9)| "incident spectrally integrated diffuse PAR"
    "net radiation of shaded soil `[W m^-2]`"
    RnSoil_diffuse::FT = 0 ;                         # | W m^-2               | (0.0, 2500)  | "net radiation of shaded soil"
    "net Short-wave radiation of sunlit soil `[W m^-2]`"
    RnSoil_direct::FT = 0;
    "net Short-wave radiation of soil (shaded + sunlit) `[W m^-2]`"
    RnSoil::FT = 0;                           # | W m^-2               | (0.0, 2500)  | "net radiation
    "net long-wave radiation of soil (shaded + sunlit) `[W m^-2]`"
    RnSoilLW::FT = 0;                           # | W m^-2               | (0.0, 2500)  | "net radiation of sunlit soil"
    # Dim of nLayers
    "net PAR of shaded leaves `[moles m^-2 s^-1]`"
    absPAR_shade::Array{FT,1} = zeros(nLayers)                  # | moles m^-2 s^-1      | (0.0, 2.5e-9)| "net PAR of shaded leaves"
    "net PAR by Cab+Car of shaded leaves `[moles m^-2 s^-1]`"
    absPAR_shadeCab::Array{FT,1} = zeros(nLayers)

    # Dimension of wavelength only:
    "Short-wave TOC outgoing radiance in observation direction `[mW m^-2 μm^-1 sr^-1]`"
    Lo::Array{FT,1} = zeros(nWl)                   # | mW m^-2 μm^-1 sr^-1  | (0.0, 2500)  | "TOC outgoing radiance in observation direction"
    "Short-wave TOC outgoing radiation `[mW m^-2 μm^-1]`"
    Eout::Array{FT,1} = zeros(nWl)                 # | mW m^-2 μm^-1        | (0.0, 2500)  | "TOC outgoing radiation
    "Short-wave Albedo in viewing direction"
    alb_obs::Array{FT,1} = zeros(nWl)              # |                      | (0.0, 1.0)  |  "albedo in viewing direction"
    "Short-wave Albedo for direct incoming radiation"
    alb_direct::Array{FT,1} = zeros(nWl)           # | -                    | (0.0, 1.0)   | "Albedo for direct incoming radiation"
    "Short-wave Albedo for diffuse incoming radiation"
    alb_diffuse::Array{FT,1} = zeros(nWl)          # | -                    | (0.0, 1.0)   | "Albedo for diffuse incoming radiation"

    # Dimension of nLevel * nWavelengths
    "upwelling diffuse short-wave radiation within canopy `[mW m^-2 μm^-1]`"
    E_up::Array{FT,2} = zeros(nWl,nLayers+1)
    "downwelling diffuse short-wave radiation within canopy `[mW m^-2 μm^-1]`"
    E_down::Array{FT,2} = zeros(nWl,nLayers+1)

    # Dimension of nLayer * nWavelengths
    "net absorbed direct radiation in each layer `[mW m^-2 μm^-1]`"
    netSW_sunlit::Array{FT,2} = zeros(nWl,nLayers)
    "net absorbed diffuse radiation in each layer `[mW m^-2 μm^-1]`"
    netSW_shade::Array{FT,2} = zeros(nWl,nLayers)
    "spectrally integrated net absorbed direct radiation in each layer `[W m^-2)]`"
    intNetSW_sunlit::Array{FT,1} = zeros(nLayers)
    "spectrally integrated net absorbed diffuse radiation in each layer `[W m^-2)]`"
    intNetSW_shade::Array{FT,1} = zeros(nLayers)
    "spectrally integrated net absorbed direct radiation in each layer `[W m^-2)]`"
    intNetLW_sunlit::Array{FT,1} = zeros(nLayers)
    "spectrally integrated net absorbed diffuse radiation in each layer `[W m^-2)]`"
    intNetLW_shade::Array{FT,1} = zeros(nLayers)


    # Dimension of nLeafInclination * nLeafAzimuth * nLayer
    "net PAR of sunlit leaves `[mol m^-2 s^-1]`"
    absPAR_sun::Array{FT,3} = zeros(nIncl,nAzi,nLayers)                    # | moles m^-2 s^-1      | (0.0, 2.5e-9)| "net PAR of sunlit leaves"
    "net PAR by Cab+Car of sunlit leaves `[mol m^-2 s^-1]`"
    absPAR_sunCab::Array{FT,3} = zeros(nIncl,nAzi,nLayers)
    "Leaf temperature (sunlit) (K)"
    T_sun3D::Array{FT,3} = zeros(nIncl,nAzi,nLayers).+285
    "Leaf temperature (sunlit) (K)"
    T_sun::Array{FT,1} = zeros(nLayers).+280
    "Fluorescence yield for sunlit leaves"
    ϕ_sun::Array{FT,3} = 0.01*ones(nIncl,nAzi,nLayers)
    "Leaf temperature (shaded) (K)"
    T_shade::Array{FT,1} = zeros(nLayers).+280
    "Fluorescence yield for shaded leaves"
    ϕ_shade::Array{FT,1} = ones(nLayers)*0.01
    "GPP of sunlit leaves `[μmol m^-2 s^-1]`"
    GPP_sun::Array{FT,3} = zeros(nIncl,nAzi,nLayers)
    "GPP of shaded leaves `[μmol m^-2 s^-1]`"
    GPP_shade::Array{FT,1} = zeros(nLayers)
    "gs of shaded leaves `[mol m^-2 s^-1]`"
    gs_shade::Array{FT,1} = zeros(nLayers)
    "gs of sunlit leaves `[mol m^-2 s^-1]`"
    gs_sun::Array{FT,3} = zeros(nIncl,nAzi,nLayers)
    "Sensible Heat flux H of shaded leaves `[W m-2]`"
    H_shade::Array{FT,1} = zeros(nLayers)
    "Sensible Heat flux H of sunlit leaves `[W m-2]`"
    H_sun::Array{FT,3} = zeros(nIncl,nAzi,nLayers)
    "Latent Heat flux LE of shaded leaves `[W m-2]`"
    LE_shade::Array{FT,1} = zeros(nLayers)
    "Latent Heat flux LE of sunlit leaves `[W m-2]`"
    LE_sun::Array{FT,3} = zeros(nIncl,nAzi,nLayers)
    "Leaf water potential of shaded leaves `[MPa]`"
    ψl_shade::Array{FT,1} = zeros(nLayers)
    "Leaf water potential of sunlit leaves `[MPa]`"
    ψl_sun::Array{FT,3} = zeros(nIncl,nAzi,nLayers)
    "Cc of shaded leaves `[µmol/mol]`"
    Cc_shade::Array{FT,1} = zeros(nLayers)
    "Cc of sunlit leaves `[µmol/mol]`"
    Cc_sun::Array{FT,3} = zeros(nIncl,nAzi,nLayers)                       
    "NPQ of shaded leaves `[-]`"
    NPQ_shade::Array{FT,1} = zeros(nLayers)
    "NPQ of sunlit leaves `[-]`"
    NPQ_sun::Array{FT,3} = zeros(nIncl,nAzi,nLayers)  

    # Fluorescence Output:
    "Hemispheric total outgoing SIF flux `[mW m^-2 μm^-1]`)"
    SIF_hemi::Array{FT,1} = zeros(nWlF)
    "Observer-direction outgoing SIF radiance  (mW m^-2 μm^-1 sr^-1))"
    SIF_obs::Array{FT,1} = zeros(nWlF)
    "Observer-direction outgoing SIF radiance, sunlit leaves  (mW m^-2 μm^-1 sr^-1))"
    SIF_obs_sunlit::Array{FT,1} = zeros(nWlF)
    "Observer-direction outgoing SIF radiance, shaded leaves  (mW m^-2 μm^-1 sr^-1))"
    SIF_obs_shaded::Array{FT,1} = zeros(nWlF)
    "Observer-direction outgoing SIF radiance, scattered   (mW m^-2 μm^-1 sr^-1))"
    SIF_obs_scattered::Array{FT,1} = zeros(nWlF)
    "Observer-direction outgoing SIF radiance, soil-reflected  (mW m^-2 μm^-1 sr^-1))"
    SIF_obs_soil::Array{FT,1} = zeros(nWlF)
    "Total SIF sum of layer sources  `[mW m^-2 μm^-1]`)"
    SIF_sum::Array{FT,1} = zeros(nWlF)

end;

"""
    struct_canopyOptProps

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct struct_canopyOptProps{FT<:AbstractFloat,nWL,nLayer,nLevel,nAzi,nIncl,nLayer_nWL, nLev_nWL,nIncl_nAzi}
    "Solar -> Diffuse backscatter weight"
    sdb::FT = 0;
    "Solar -> Diffuse forward scatter weight"
    sdf::FT = 0;
    "Diffuse -> Directional backscatter weight"
    dob::FT = 0;
    "Diffuse -> Directional forward scatter weight"
    dof::FT = 0;
    "Diffuse -> Diffuse backscatter weight"
    ddb::FT = 0;
    "Diffuse -> Diffuse forward scatter weight"
    ddf::FT = 0;
    "Solar beam extinction coefficient weight"
    ks::FT = 0;
    "Outgoing beam extinction coefficient weight"
    ko::FT = 0;
    ""
    bf::FT = 0;
    "weight of specular2directional backscatter coefficient"
    sob::FT = 0;
    "weight of specular2directional forward coefficient"
    sof::FT = 0;

    # now multi dimensional arrays:
    "conversion factor fs to compute irradiance on inclined leaf"
    fs::MArray{Tuple{nIncl, nAzi}, FT,2,nIncl_nAzi}          = MArray{Tuple{nIncl, nAzi}, FT}(undef);
    "abs(fs)"
    absfs::MArray{Tuple{nIncl, nAzi}, FT,2,nIncl_nAzi}     = MArray{Tuple{nIncl, nAzi}, FT}(undef);
    "abs(fs*fo)"
    absfsfo::MArray{Tuple{nIncl, nAzi}, FT,2,nIncl_nAzi}    = MArray{Tuple{nIncl, nAzi}, FT}(undef);
    "fs*fo"
    fsfo::MArray{Tuple{nIncl, nAzi}, FT,2,nIncl_nAzi}          = MArray{Tuple{nIncl, nAzi}, FT}(undef);
    "conversion factor fo for angle towards observer (not sun like fs)"
    fo::MArray{Tuple{nIncl, nAzi}, FT,2,nIncl_nAzi}          = MArray{Tuple{nIncl, nAzi}, FT}(undef);
    "Cosine of leaf azimuths"
    cosΘ_l ::MArray{Tuple{nIncl, nAzi}, FT,2,nIncl_nAzi}          = MArray{Tuple{nIncl, nAzi}, FT}(undef);
    "cos of leaf azimuth sqared"
    cos2Θ_l ::MArray{Tuple{nIncl, nAzi}, FT,2,nIncl_nAzi}           = MArray{Tuple{nIncl, nAzi}, FT}(undef);
    "Probability of directly viewing a leaf in solar direction"
    Ps::MArray{Tuple{nLevel}, FT,1,nLevel}           = MArray{Tuple{nLevel}, FT}(undef);
    "Probability of directly viewing a leaf in viewing direction"
    Po::MArray{Tuple{nLevel}, FT,1,nLevel}           = MArray{Tuple{nLevel}, FT}(undef);
    "Bi-directional probability of directly viewing a leaf (solar->canopy->viewing)"
    Pso::MArray{Tuple{nLevel}, FT,1,nLevel}            = MArray{Tuple{nLevel}, FT}(undef);

    # The following also depend on leaf reflectance and transmission. Might go into a separate strcuture so that we can have it separately for thermal, SW and SIF?
    "diffuse     backscatter scattering coefficient for diffuse  incidence"
    sigb::MArray{Tuple{nWL, nLayer}, FT,2,nLayer_nWL}             = MArray{Tuple{nWL, nLayer}, FT}(undef);
    "diffuse     forward     scattering coefficient for diffuse  incidence"
    sigf::MArray{Tuple{nWL, nLayer}, FT,2,nLayer_nWL}             = MArray{Tuple{nWL, nLayer}, FT}(undef);
    "diffuse     backscatter scattering coefficient for specular incidence"
    sb::MArray{Tuple{nWL, nLayer}, FT,2,nLayer_nWL}               = MArray{Tuple{nWL, nLayer}, FT}(undef);
    "diffuse     forward     scattering coefficient for specular incidence"
    sf::MArray{Tuple{nWL, nLayer}, FT,2,nLayer_nWL}                 = MArray{Tuple{nWL, nLayer}, FT}(undef);
    "directional backscatter scattering coefficient for diffuse  incidence"
    vb::MArray{Tuple{nWL, nLayer}, FT,2,nLayer_nWL}                = MArray{Tuple{nWL, nLayer}, FT}(undef);
    "directional forward     scattering coefficient for diffuse  incidence"
    vf::MArray{Tuple{nWL, nLayer}, FT,2,nLayer_nWL}                = MArray{Tuple{nWL, nLayer}, FT}(undef);
    "bidirectional scattering coefficent (directional-directional)"
    w::MArray{Tuple{nWL, nLayer}, FT,2,nLayer_nWL}                 = MArray{Tuple{nWL, nLayer}, FT}(undef);
    "attenuation"
    a::MArray{Tuple{nWL, nLayer}, FT,2,nLayer_nWL}                  = MArray{Tuple{nWL, nLayer}, FT}(undef);
    "Effective layer transmittance (direct->diffuse)"
    Xsd::MArray{Tuple{nWL, nLayer}, FT,2,nLayer_nWL}             = MArray{Tuple{nWL, nLayer}, FT}(undef);
    "Effective layer transmittance (diffuse->diffuse)"
    Xdd::MArray{Tuple{nWL, nLayer}, FT,2,nLayer_nWL}              = MArray{Tuple{nWL, nLayer}, FT}(undef);
    "Effective layer reflectance (direct->diffuse)"
    R_sd::MArray{Tuple{nWL, nLevel}, FT,2,nLev_nWL}         = MArray{Tuple{nWL, nLevel}, FT}(undef);
    "Effective layer reflectance (diffuse->diffuse)"
    R_dd::MArray{Tuple{nWL, nLevel}, FT,2,nLev_nWL}        = MArray{Tuple{nWL, nLevel}, FT}(undef);

    "Solar direct radiation per layer)"
    Es_::MArray{Tuple{nWL, nLevel}, FT,2,nLev_nWL}        = MArray{Tuple{nWL, nLevel}, FT}(undef);
end;

function create_canopyOpt(; FType=FT, nWL::Int=114, nLayers::Int=20, nAzi::Int=36, nIncl::Int=9)
    return struct_canopyOptProps{FType,nWL,nLayers,nLayers+1,nAzi,nIncl,nLayers*nWL, (nLayers+1)*nWL,nIncl*nAzi }();
end

"""
    struct_canopy

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct struct_canopy{FT<:AbstractFloat, n_layer, lai}
    "number of canopy layers" # This needs to come globally though, can lead to errors right now as defined elsewhere as well!
    nlayers::Int     = n_layer    # | -          | (2, 60)      | "number of canopy layers"
    "Leaf Area Index"
    LAI::FT            = lai   # | -          | (0.0, 9.0)   | "Leaf Area Index"
    "Clumping factor"
    Ω::FT             = 1.0   # | -          | (0.0, 1.0)   | "clumping factor"
    "Structure factor a"
    clump_a::FT       = 1.0   # | -          | (0.0, 1.0)   | "structure factor a"
    "Structure factor b"
    clump_b::FT       = 0.0   # | -          | (0.0, 1.0)   | "structure factor b (angular variation)"       
    "Leaf width"
    leafwidth::FT      = 0.1   # | m          | (0.0, 1.0)   | "Leaf width"
    "Vegetation height"
    hc::FT             = 2.0   # | m          | (0.0, 70.0)  | "Vegetation height"
    "Leaf Inclination"
    LIDFa::FT          = 0. # | -          | (-1.0, 1.0)  | "Leaf Inclination"
    "Variation in leaf inclination"
    LIDFb::FT          = 0. # | -          | (-1.0, 1.0)  | "Variation in leaf inclination"
    "HotSpot parameter (still need to check!)"
    hot::FT            = 0.05  # | -          | (0, 1.0)     | "HotSpot parameter (still need to check!)"
    # tree/leaf traits
    "Canopy height (m)"
    height::FT      = 20.;                            # tree height (m)
    "Canopy roughness (m)"
    z0m::FT         = 1.;                          # tree roughness (m)
    z0h::FT         = -999.;                          # tree roughness (m) - TODO should be changed later
    "Canopy displacement height (m)"
    d::FT           = -999.;                          # tree displacement height (m)
    "m/sqrt(s) turbulent transfer coefficient"
    Cd::FT          = 0.01;  
    

    # Some more derived parameters:
    lazitab::Array{FT} = collect(5.0:10.0:355.0)
    # This is changed afterwards, ignore here.
    lidf::Array{FT} = similar(collect(5.0:10.0:355.0))
    xl::Array{FT} = collect(0.0:-1.0/nlayers:-1.0)
    dx::FT             = 1.0/nlayers
end


function loadOpti(swl::Array; file=file_Opti)
    # Read in all optical data:
    FT = eltype(swl)
    opti = matread(file)["optipar"]
    nr_     =  opti["nr"]
    Km_     =  opti["Kdm"]
    Kab_    =  opti["Kab"]
    Kant_   =  opti["Kant"]
    Kcar_   =  opti["Kca"]
    Kw_     =  opti["Kw"]
    KBrown_ =  opti["Ks"]
    phi_    =  opti["phi"]
    KcaV_   =  opti["KcaV"]
    KcaZ_   =  opti["KcaZ"]
    lambda_ =  opti["wl"]
    #@show length(swl)
    nr = MArray{Tuple{length(swl)-1}, FT}(undef)
    Km =MArray{Tuple{length(swl)-1}, FT}(undef)
    Kab = MArray{Tuple{length(swl)-1}, FT}(undef)
    Kant = MArray{Tuple{length(swl)-1}, FT}(undef)
    Kcar = MArray{Tuple{length(swl)-1}, FT}(undef)
    Kw =MArray{Tuple{length(swl)-1}, FT}(undef)
    KBrown = MArray{Tuple{length(swl)-1}, FT}(undef)
    phi = MArray{Tuple{length(swl)-1}, FT}(undef)
    KcaV = MArray{Tuple{length(swl)-1}, FT}(undef)
    KcaZ = MArray{Tuple{length(swl)-1}, FT}(undef)
    lambda = MArray{Tuple{length(swl)-1}, FT}(undef)
    kChlrel = MArray{Tuple{length(swl)-1}, FT}(undef)
    println("Reading Optical Parameters from ", swl[1], " to ", swl[end], " length: ", length(swl))
    @inbounds for i in 1:length(swl)-1
        wo = findall((lambda_.>=swl[i]).&(lambda_.<swl[i+1]) )
        if length(wo)==0
            println("Warning, some wavelengths out of bounds ", swl[i])
        end
        #@show typeof(mean(nr_[wo]))
        nr[i]     =  mean(nr_[wo])
        Km[i]     =  mean(Km_[wo])
        Kab[i]    =  mean(Kab_[wo])
        Kant[i]   =  mean(Kant_[wo])
        Kcar[i]   =  mean(Kcar_[wo])
        Kw[i]     =  mean(Kw_[wo])
        KBrown[i] =  mean(KBrown_[wo])
        phi[i]    =  mean(phi_[wo])
        KcaV[i]   =  mean(KcaV_[wo])
        KcaZ[i]   =  mean(KcaZ_[wo])
        lambda[i] =  mean(lambda_[wo])
    end
    return nr, Km, Kab, Kant, Kcar, Kw, KBrown, phi, KcaV, KcaZ, lambda
end


function loadSun(swl::Array; file=file_Sun)
    FT = eltype(swl)
    # Read in all optical data:
    suni  = matread(file)["sun"]
    wl    =  suni["wl"]
    Edir  =  suni["Edirect"]
    Ediff =  suni["Ediffuse"]

    wl_    = MArray{Tuple{length(swl)-1}, Float32}(undef)
    Edir_  = MArray{Tuple{length(swl)-1}, Float32}(undef)
    Ediff_ = MArray{Tuple{length(swl)-1}, Float32}(undef)
    #println("Reading Optical Parameters from ", swl[1], " to ", swl[end], " length: ", length(swl))
    for i in 1:length(swl)-1
        wo = findall((wl.>=swl[i]).&(wl.<swl[i+1]) )
        if length(wo)==0
            println("Warning, some wavelengths out of bounds ", swl[i])
        end
        wl_[i]   =  mean(wl[wo])
        Edir_[i]  =  mean(Edir[wo])
        Ediff_[i]  = mean(Ediff[wo])
    end
    return wl_, Edir_, Ediff_
end
