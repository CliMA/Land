using Parameters
using DocStringExtensions
using StaticArrays

# Fixed path right now here
file_Opti = joinpath(dirname(pathof(CanopyRTMod)), "Optipar2017_ProspectD.mat")
file_Sun = joinpath(dirname(pathof(CanopyRTMod)), "sun.mat")

# Struct for observation and solar angles
"""
    Struct for observation and solar angles

# Fields
$(DocStringExtensions.FIELDS)
"""
@with_kw mutable struct struct_angles{FT<:Number}
    "Solar Zenith Angle in degrees"
    tts::FT=30
    "Viewing Zenith Angle in degrees"
    tto::FT=0
    "relative azimuth in degrees"
    psi::FT=0
end

# Struct for leaf optical properties
struct optipar{FT<:Number}
    nr::Array{FT}
    Km::Array{FT}
    Kab::Array{FT}
    Kant::Array{FT}
    Kcar::Array{FT}
    Kw::Array{FT}
    KBrown::Array{FT}
    phi::Array{FT}
    KcaV::Array{FT}
    KcaZ::Array{FT}
    lambda::Array{FT}
end

"""
    incomingRadiation

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct incomingRadiation{FT<:Number}
    "Wavelength (nm)"
    wl::Array{FT,1}
    "Direct incoming radiation (mW m^-2 μm^-1)"
    E_direct::Array{FT,1}
    " Diffuse incoming radiation (mW m^-2 μm^-1)"
    E_diffuse::Array{FT,1}

end

"""
    leafbio

# Fields
$(DocStringExtensions.FIELDS)
"""
@with_kw mutable struct leafbio{FT<:Number,nWl,nWle, nWlf}
    "Leaf structure parameter"
    N::FT    = 1.4       # | -          | (1.0, 3.0)  | "Leaf structure parameter"
    "Chlorophyll a+b content"
    Cab::FT  = 40.0      # | μg cm^-2   | (0.0, 110)  | "Chlorophyll a+b content"
    "Carotenoid content"
    Car::FT  = 10.0      # | μg cm^-2   | (0.0, 40.0) | "Carotenoid content"
    "Anthocynanin content"
    Ant::FT  = 0.0       # | μg cm^-2   | (0.0, 40.0) | "Anthocynanin content"
    "Senescent material fraction"
    Cs::FT   = 0.0       # | -          | (0.0, 1.0)  | "Senescent material fraction"
    "Equivalent water thickness"
    Cw::FT   = 0.009     # | cm         | (0.0, 0.05) | "Equivalent water thickness"
    "Dry matter content (dry leaf mass per unit area)"
    Cm::FT   = 0.012      # | g cm^-2    | (0.0, 0.2)  | "Dry matter content (dry leaf mass per unit area)"
    "Fractionation between Zeaxanthin and Violaxanthin in Car (1=all Zeaxanthin)"
    Cx::FT   = 0.0       # | -          | (0.0, 1.0)  | "Fractionation between Zeaxanthin and Violaxanthin in Car (1=all Zeaxanthin)"
    "Broadband thermal reflectance"
    ρ_LW::FT = 0.01      # | -          | (0.0, 1.0)  | "Broadband thermal reflectance"
    "Broadband thermal transmission"
    τ_LW::FT = 0.01      # | -          | (0.0, 1.0)  | "Broadband thermal transmission"
    "Leaf fluorescence efficiency"
    fqe::FT = 0.01       # | -          | (0.0, 1.0)  | "Leaf fluorescence efficiency"
    "shortwave leaf reflectance"
    ρ_SW::Array{FT,1} = zeros(nWl)    # | -          | (0.0, 1.0)  | "shortwave reflectance"
    "shortwave leaf transmission"
    τ_SW::Array{FT,1} = zeros(nWl)    # | -          | (0.0, 1.0)  | "shortwave transmission"
    "relative absorbtion by Chlorophyll+Car"
    kChlrel::Array{FT,1} = zeros(nWl) # | -          | (0.0, 1.0)  | "relative absorbtion by Chlorophyll"
    "relative absorbtion by Chlorophyll"
    kChlrel_old::Array{FT,1} = zeros(nWl) # | -          | (0.0, 1.0)  | "relative absorbtion by Chlorophyll"
    "Fluorescence excitation matrix backwards"
    Mb::Array{FT,2}= zeros(nWle,nWlf)      # | -          | (0.0, 1.0)  | "Fluorescence excitation matrix backwards"
    "Fluorescence excitation matrix forwards"
    Mf::Array{FT,2}= zeros(nWle,nWlf)      # | -          | (0.0, 1.0)  | "Fluorescence excitation matrix forwards"

end

"""
    struct_soil

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct struct_soil{FT<:Number}
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
@with_kw mutable struct struct_canopyRadiation{FT<:AbstractFloat,nWl,nWlF,nIncl,nAzi,nLayers}
    # Scalars
    "integrated TOC outgoing flux (W m^-2)"
    intEout::FT = 0;                       # | W m^-2               | (0.0, 2500)  | "integrated TOC outgoing flux"
    "incident spectrally integrated total PAR (moles m^-2 s^-1)"
    incomingPAR::FT = 0;                           # | moles m^-2 s^-1      | (0.0, 2.5e-9)| "incident spectrally integrated total PAR"
    "incident spectrally integrated direct PAR (moles m^-2 s^-1)"
    incomingPAR_direct::FT = 0;                    # | moles m^-2 s^-1      | (0.0, 2.5e-9)| "incident spectrally integrated direct PAR"
    "incident spectrally integrated diffuse PAR (moles m^-2 s^-1)"
    incomingPAR_diffuse::FT = 0;                  # | moles m^-2 s^-1      | (0.0, 2.5e-9)| "incident spectrally integrated diffuse PAR"
    "net radiation of shaded soil (W m^-2)"
    RnSoil_diffuse::FT = 0 ;                         # | W m^-2               | (0.0, 2500)  | "net radiation of shaded soil"
    "net Short-wave radiation of sunlit soil (W m^-2)"
    RnSoil_direct::FT = 0;
    "net Short-wave radiation of soil (shaded + sunlit) (W m^-2)"
    RnSoil::FT = 0;                           # | W m^-2               | (0.0, 2500)  | "net radiation
    "net long-wave radiation of soil (shaded + sunlit) (W m^-2)"
    RnSoilLW::FT = 0;                           # | W m^-2               | (0.0, 2500)  | "net radiation of sunlit soil"
    # Dim of nLayers
    "net PAR of shaded leaves (moles m^-2 s^-1)"
    absPAR_shade::Array{FT,1} = zeros(nLayers)                  # | moles m^-2 s^-1      | (0.0, 2.5e-9)| "net PAR of shaded leaves"
    "net PAR by Cab+Car of shaded leaves (moles m^-2 s^-1)"
    absPAR_shadeCab::Array{FT,1} = zeros(nLayers)

    # Dimension of wavelength only:
    "Short-wave TOC outgoing radiance in observation direction (mW m^-2 μm^-1 sr^-1)"
    Lo::Array{FT,1} = zeros(nWl)                   # | mW m^-2 μm^-1 sr^-1  | (0.0, 2500)  | "TOC outgoing radiance in observation direction"
    "Short-wave TOC outgoing radiation (mW m^-2 μm^-1)"
    Eout::Array{FT,1} = zeros(nWl)                 # | mW m^-2 μm^-1        | (0.0, 2500)  | "TOC outgoing radiation
    "Short-wave Albedo in viewing direction"
    alb_obs::Array{FT,1} = zeros(nWl)              # |                      | (0.0, 1.0)  |  "albedo in viewing direction"
    "Short-wave Albedo for direct incoming radiation"
    alb_direct::Array{FT,1} = zeros(nWl)           # | -                    | (0.0, 1.0)   | "Albedo for direct incoming radiation"
    "Short-wave Albedo for diffuse incoming radiation"
    alb_diffuse::Array{FT,1} = zeros(nWl)          # | -                    | (0.0, 1.0)   | "Albedo for diffuse incoming radiation"

    # Dimension of nLayer+1 * nWavelengths
    "upwelling diffuse short-wave radiation within canopy (mW m^-2 μm^-1)"
    E_up::Array{FT,2} = zeros(nWl,nLayers+1)
    "downwelling diffuse short-wave radiation within canopy (mW m^-2 μm^-1)"
    E_down::Array{FT,2} = zeros(nWl,nLayers+1)

    # Dimension of nLayer * nWavelengths
    "net absorbed direct radiation in each layer (mW m^-2 μm^-1)"
    netSW_sunlit::Array{FT,2} = zeros(nWl,nLayers)
    "net absorbed diffuse radiation in each layer (mW m^-2 μm^-1)"
    netSW_shade::Array{FT,2} = zeros(nWl,nLayers)
    "spectrally integrated net absorbed direct radiation in each layer (W m^-2)"
    intNetSW_sunlit::Array{FT,1} = zeros(nLayers)
    "spectrally integrated net absorbed diffuse radiation in each layer (W m^-2)"
    intNetSW_shade::Array{FT,1} = zeros(nLayers)
    "spectrally integrated net absorbed direct radiation in each layer (W m^-2)"
    intNetLW_sunlit::Array{FT,1} = zeros(nLayers)
    "spectrally integrated net absorbed diffuse radiation in each layer (W m^-2)"
    intNetLW_shade::Array{FT,1} = zeros(nLayers)


    # Dimension of nLeafInclination * nLeafAzimuth * nLayer
    "net PAR of sunlit leaves moles m^-2 s^-1"
    absPAR_sun::Array{FT,3} = zeros(nIncl,nAzi,nLayers)                    # | moles m^-2 s^-1      | (0.0, 2.5e-9)| "net PAR of sunlit leaves"
    "net PAR by Cab+Car of sunlit leaves moles m^-2 s^-1"
    absPAR_sunCab::Array{FT,3} = zeros(nIncl,nAzi,nLayers)
    "Leaf temperature (sunlit) (K)"
    T_sun3D::Array{FT,3} = zeros(nIncl,nAzi,nLayers).+285
    "Leaf temperature (sunlit) (K)"
    T_sun::Array{FT,1} = zeros(nLayers).+280
    "Fluorescence yield for sunlit leaves"
    φ_sun::Array{FT,3} = ones(nIncl,nAzi,nLayers)
    "Leaf temperature (shaded) (K)"
    T_shade::Array{FT,1} = zeros(nLayers).+280
    "Fluorescence yield for shaded leaves"
    φ_shade::Array{FT,1} = ones(nLayers)

    # Fluorescence Output:
    "Hemispheric total outgoing SIF flux (mW m^-2 μm^-1))"
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
    "Total SIF sum of layer sources  (mW m^-2 μm^-1))"
    SIF_sum::Array{FT,1} = zeros(nWlF)

end

"""
    struct_canopyRadiation

# Fields
$(DocStringExtensions.FIELDS)
"""
@with_kw mutable struct struct_canopyOptProps{FT<:AbstractFloat,nWL,nLayer,nAzi,nIncl}
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
    "per leaf angles"
    fs::Array{FT,2} = zeros(nIncl, nAzi);
    "per leaf angles"
    absfs::Array{FT,2} = zeros(nIncl, nAzi);
    "abs(fs*fo)"
    absfsfo::Array{FT,2} = zeros(nIncl, nAzi);
    "fs*fo"
    fsfo::Array{FT,2} = zeros(nIncl, nAzi);
    "per leaf angles"
    fo::Array{FT,2} = zeros(nIncl, nAzi);
    "Cosine of leaf azimuths"
    cosΘ_l ::Array{FT,2} = zeros(nIncl, nAzi);
    "cos of leaf azimuth sqared"
    cos2Θ_l ::Array{FT,2} = zeros(nIncl, nAzi);
    "Probability of directly viewing a leaf in solar direction"
    Ps::Array{FT,1} = zeros(nLayer+1)
    "Probability of directly viewing a leaf in viewing direction"
    Po::Array{FT,1} = zeros(nLayer+1)
    "Bi-directional probability of directly viewing a leaf (solar->canopy->viewing)"
    Pso::Array{FT,1} = zeros(nLayer+1)

    # The following also depend on leaf reflectance and transmission. Might go into a separate strcuture so that we can have it separately for thermal, SW and SIF?
    "diffuse     backscatter scattering coefficient for diffuse  incidence"
    sigb::Array{FT,2} = zeros(nWL, nLayer)
    "diffuse     forward     scattering coefficient for diffuse  incidence"
    sigf::Array{FT,2} = zeros(nWL, nLayer)
    "diffuse     backscatter scattering coefficient for specular incidence"
    sb::Array{FT,2} = zeros(nWL, nLayer)
    "diffuse     forward     scattering coefficient for specular incidence"
    sf::Array{FT,2} = zeros(nWL, nLayer)
    "directional backscatter scattering coefficient for diffuse  incidence"
    vb::Array{FT,2} = zeros(nWL, nLayer)
    "directional forward     scattering coefficient for diffuse  incidence"
    vf::Array{FT,2} = zeros(nWL, nLayer)
    "bidirectional scattering coefficent (directional-directional)"
    w::Array{FT,2} = zeros(nWL, nLayer)
    "attenuation"
    a::Array{FT,2} = zeros(nWL, nLayer)
    "Effective layer transmittance (direct->diffuse)"
    Xsd::Array{FT,2} = zeros(nWL, nLayer)
    "Effective layer transmittance (diffuse->diffuse)"
    Xdd::Array{FT,2} = zeros(nWL, nLayer)
    "Effective layer reflectance (direct->diffuse)"
    R_sd::Array{FT,2} = zeros(nWL, nLayer+1)
    "Effective layer reflectance (diffuse->diffuse)"
    R_dd::Array{FT,2} = zeros(nWL, nLayer+1)

    "Solar direct radiation per layer)"
    Es_::Array{FT,2} = zeros(nWL, nLayer+1)

end


"""
    struct_canopy

# Fields
$(DocStringExtensions.FIELDS)
"""
@with_kw mutable struct struct_canopy{FT<:Number}
    "number of canopy layers" # This needs to come globally though, can lead to errors right now as defined elsewhere as well!
    nlayers::Int     = 20    # | -          | (2, 60)      | "number of canopy layers"
    "Leaf Area Index"
    LAI::FT            = 3.0   # | -          | (0.0, 9.0)   | "Leaf Area Index"
    "Clumping factor"
    Ω::FT             = 1.0   # | -          | (0.0, 1.0)   | "clumping factor"
    "Leaf width"
    leafwidth::FT      = 0.1   # | m          | (0.0, 1.0)   | "Leaf width"
    "Vegetation height"
    hc::FT             = 2.0   # | m          | (0.0, 70.0)  | "Vegetation height"
    "Leaf Inclination"
    LIDFa::FT          = -0.35 # | -          | (-1.0, 1.0)  | "Leaf Inclination"
    "Variation in leaf inclination"
    LIDFb::FT          = -0.15 # | -          | (-1.0, 1.0)  | "Variation in leaf inclination"
    "HotSpot parameter (still need to check!)"
    hot::FT            = 0.05  # | -          | (0, 1.0)     | "HotSpot parameter (still need to check!)"
    "Leaf distribution type (2=campbell, 1=ladgen)"
    TypeLidf::FT       = 1     # | -          | (-1.0, 1.0)  | "Leaf distribution type (2=campbell, 1=ladgen)"

    # Some more derived parameters:
    lazitab::Array{FT} = collect(5.0:10.0:355.0)
    # This is changed afterwards, ignore here.
    lidf::Array{FT} = similar(collect(5.0:10.0:355.0))
    xl::Array{FT} = collect(0.0:-1.0/nlayers:-1.0)
    dx::FT             = 1.0/nlayers
end


function loadOpti(swl::AbstractArray; file=file_Opti)
    # Read in all optical data:
    FT = typeof(swl)
    opti = matread(file_Opti)["optipar"]
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

    nr = FT(undef,length(swl)-1)
    Km = FT(undef,length(swl)-1)
    Kab = FT(undef,length(swl)-1)
    Kant = FT(undef,length(swl)-1)
    Kcar = FT(undef,length(swl)-1)
    Kw = FT(undef,length(swl)-1)
    KBrown = FT(undef,length(swl)-1)
    phi = FT(undef,length(swl)-1)
    KcaV = FT(undef,length(swl)-1)
    KcaZ = FT(undef,length(swl)-1)
    lambda = FT(undef,length(swl)-1)
    kChlrel = FT(undef,length(swl)-1)
    println("Reading Optical Parameters from ", swl[1], " to ", swl[end], " length: ", length(swl))
    for i in 1:length(swl)-1
        wo = findall((lambda_.>=swl[i]).&(lambda_.<swl[i+1]) )
        if length(wo)==0
            println("Warning, some wavelengths out of bounds ", swl[i])
        end
        nr[i]   =  mean(nr_[wo])
        Km[i]  =  mean(Km_[wo])
        Kab[i]  =  mean(Kab_[wo])
        Kant[i] =  mean(Kant_[wo])
        Kcar[i]  =  mean(Kcar_[wo])
        Kw[i]   =  mean(Kw_[wo])
        KBrown[i] = mean(KBrown_[wo])
        phi[i]  =  mean(phi_[wo])
        KcaV[i] = mean(KcaV_[wo])
        KcaZ[i] = mean(KcaZ_[wo])
        lambda[i] = mean(lambda_[wo])
    end
    return nr, Km, Kab, Kant, Kcar, Kw, KBrown, phi, KcaV, KcaZ, lambda
end


function loadSun(swl::AbstractArray, file=file_Sun)
    FT = typeof(swl)
    # Read in all optical data:
    suni = matread(file_Sun)["sun"]
    wl   =  suni["wl"]
    Edir =  suni["Edirect"]
    Ediff =  suni["Ediffuse"]

    wl_ = FT(undef,length(swl)-1)
    Edir_ = FT(undef,length(swl)-1)
    Ediff_ = FT(undef,length(swl)-1)
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
