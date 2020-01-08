using Parameters
using DocStringExtensions

# Fixed path right now here
file_Opti = joinpath(dirname(pathof(FluspectMod)), "Optipar2017_ProspectD.mat")
file_Sun = joinpath(dirname(pathof(FluspectMod)), "sun.mat")

# Struct for observation and solar angles
"""
    Struct for observation and solar anglesn

# Fields
$(DocStringExtensions.FIELDS)
"""
@with_kw mutable struct struct_angles{FT<:Number}
    "Solar Zenith Angle in degrees"
    tts::FT=45
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
    struct_canopyRadiation

# Fields
$(DocStringExtensions.FIELDS)
"""
@with_kw mutable struct leafbio{FT<:Number}
    "Leaf structure parameter"
    N::FT    = 1.5       # | -          | (1.0, 3.0)  | "Leaf structure parameter"
    "Chlorophyll a+b content"
    Cab::FT  = 40.0      # | μg cm^-2   | (0.0, 110)  | "Chlorophyll a+b content"
    "Carotenoid content"
    Car::FT  = 10.0      # | μg cm^-2   | (0.0, 40.0) | "Carotenoid content"
    "Anthocynanin content"
    Ant::FT  = 8.0       # | μg cm^-2   | (0.0, 40.0) | "Anthocynanin content"
    "Senescent material fraction"
    Cs::FT   = 0.0       # | -          | (0.0, 1.0)  | "Senescent material fraction"
    "Equivalent water thickness"
    Cw::FT   = 0.015     # | cm         | (0.0, 0.05) | "Equivalent water thickness"
    "Dry matter content (dry leaf mass per unit area)"
    Cm::FT   = 0.01      # | g cm^-2    | (0.0, 0.2)  | "Dry matter content (dry leaf mass per unit area)"
    "Fractionation between Zeaxanthin and Violaxanthin in Car (1=all Zeaxanthin)"
    Cx::FT   = 0.0       # | -          | (0.0, 1.0)  | "Fractionation between Zeaxanthin and Violaxanthin in Car (1=all Zeaxanthin)"
    "Broadband thermal reflectance"
    ρ_LW::FT = 0.01      # | -          | (0.0, 1.0)  | "Broadband thermal reflectance"
    "Broadband thermal transmission"
    τ_LW::FT = 0.01      # | -          | (0.0, 1.0)  | "Broadband thermal transmission"
    "Leaf fluorescence efficiency"
    fqe::FT = 0.01       # | -          | (0.0, 1.0)  | "Leaf fluorescence efficiency"
    "shortwave leaf reflectance"
    ρ_SW::Array{FT,1}    # | -          | (0.0, 1.0)  | "shortwave reflectance"
    "shortwave leaf transmission"
    τ_SW::Array{FT,1}    # | -          | (0.0, 1.0)  | "shortwave transmission"
    "relative absorbtion by Chlorophyll"
    kChlrel::Array{FT,1} # | -          | (0.0, 1.0)  | "relative absorbtion by Chlorophyll"
    "Fluorescence excitation matrix backwards"
    Mb::Array{FT,2}      # | -          | (0.0, 1.0)  | "Fluorescence excitation matrix backwards"
    "Fluorescence excitation matrix forwards"
    Mf::Array{FT,2}      # | -          | (0.0, 1.0)  | "Fluorescence excitation matrix forwards"

end

"""
    struct_canopyRadiation

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct struct_soil{FT<:Number}
    "shortwave albedo"
    albedo_SW::Array{FT,1}    # | -          | (0.0, 1.0)  | "shortwave albedo"
    "longwave albedo"
    albedo_LW::Array{FT,1}    # | -          | (0.0, 1.0)  | "longwave albedo"
end

"""
    struct_canopyRadiation

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct struct_canopyRadiation{FT<:Number}
    # Scalars
    "integrated TOC outgoing flux"
    intEout::FT                       # | W m^-2               | (0.0, 2500)  | "integrated TOC outgoing flux"
    "incident spectrally integrated total PAR"
    PAR::FT                           # | moles m^-2 s^-1      | (0.0, 2.5e-9)| "incident spectrally integrated total PAR"
    "incident spectrally integrated direct PAR"
    PAR_direct::FT                    # | moles m^-2 s^-1      | (0.0, 2.5e-9)| "incident spectrally integrated direct PAR"
    "incident spectrally integrated diffuse PAR"
    PAR_diffuse::FT                   # | moles m^-2 s^-1      | (0.0, 2.5e-9)| "incident spectrally integrated diffuse PAR"
    "net radiation of shaded soil"
    Rnhs::FT                          # | W m^-2               | (0.0, 2500)  | "net radiation of shaded soil"
    "net radiation of sunlit soil"
    Rnus::FT                          # | W m^-2               | (0.0, 2500)  | "net radiation of sunlit soil"
    "extinction cofficient in the solar direction"
    k::FT
    "extinction cofficient in the viewing direction"
    K::FT

    # Dim of nLayers+1
    "gap fraction in the solar direction"
    Ps::Array{FT,1}
    "gap fraction in the viewing direction"
    Po::Array{FT,1}
    "bi-derectional gap fraction (solar->canopy->viewing)"
    Pso::Array{FT,1}

    # Dim of nLayers
    "net radiation of shaded leaves"
    Rnhc::Array{FT,1}                 # | W m^-2               | (0.0, 2500)  | "net radiation of shaded leaves"
    "net PAR of shaded leaves"
    Pnh::Array{FT,1}                  # | moles m^-2 s^-1      | (0.0, 2.5e-9)| "net PAR of shaded leaves"
    "net PAR absorbed by Cab of shaded leaves"
    Pnh_Cab::Array{FT,1}              # | moles m^-2 s^-1      | (0.0, 2.5e-9)| "net PAR absorbed by Cab of shaded leaves"
    "net PAR absorbed by shaded leaves"
    Rnh_PAR::Array{FT,1}              # | W m^-2               | (0.0, 2500)  | "net PAR absorbed by shaded leaves

    # Dimension of wavelength only:
    "TOC outgoing radiance in observation direction"
    Lo::Array{FT,1}                   # | mW m^-2 μm^-1 sr^-1  | (0.0, 2500)  | "TOC outgoing radiance in observation direction"
    "TOC outgoing radiation"
    Eout::Array{FT,1}                 # | mW m^-2 μm^-1        | (0.0, 2500)  | "TOC outgoing radiation"
    "incident direct radiation at top of canopy"
    inc_SW_direct::Array{FT,1}        # | mW m^-2 μm^-1        | (0.0, 2500)  | "incident direct radiation at top of canopy"
    "incident diffuse radiation at top of canopy"
    inc_SW_diffuse::Array{FT,1}       # | mW m^-2 μm^-1        | (0.0, 2500)  | "incident diffuse radiation at top of canopy"
    "Albedo for direct incoming radiation"
    alb_direct::Array{FT,1}           # | -                    | (0.0, 1.0)   | "Albedo for direct incoming radiation"
    "Albedo for diffuse incoming radiation"
    alb_diffuse::Array{FT,1}          # | -                    | (0.0, 1.0)   | "Albedo for diffuse incoming radiation"


    # Dimension of nLayer+1 * nWavelengths
    "upwelling diffuse radiation within canopy"
    E_up::Array{FT,2}                 # | mW m^-2 μm^-1        | (0.0, 2500)  | "upwelling diffuse radiation within canopy"
    "downwelling diffuse radiation within canopy"
    E_down::Array{FT,2}               # | mW m^-2 μm^-1        | (0.0, 2500)  | "downwelling diffuse radiation within canopy"

    # Dimension of nLeafInclination * nLeafAzimuth * nLayer
    "net radiation of sunlit leaves"
    Rnuc::Array{FT,3}                 # | W m^-2               | (0.0, 2500)  | "net radiation of sunlit leaves"
    "net PAR of sunlit leaves"
    Pnu::Array{FT,3}                  # | moles m^-2 s^-1      | (0.0, 2.5e-9)| "net PAR of sunlit leaves"
    "net PAR absorbed by Cab of sunlit leaves"
    Pnu_Cab::Array{FT,3}              # | moles m^-2 s^-1      | (0.0, 2.5e-9)| "net PAR absorbed by Cab of sunlit leaves"
    "net PAR absorbed by sunlit leaves"
    Rnu_PAR::Array{FT,3}              # | W m^-2               | (0.0, 2500)  | "net PAR absorbed by sunlit leaves"
end

"""
    struct_canopy

# Fields
$(DocStringExtensions.FIELDS)
"""
@with_kw mutable struct struct_canopy{FT<:Number}
    "number of canopy layers"
    nlayers::Int64     = 20    # | -          | (2, 60)      | "number of canopy layers"
    "Leaf Area Index"
    LAI::FT            = 3.0   # | -          | (0.0, 9.0)   | "Leaf Area Index"
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
    xl::Array{FT}      = collect(0.0:-1.0/nlayers:-1.0)
    dx::FT             = 1.0/nlayers
end

function getRadStruct(nWL::Int,nLayer::Int,nAzi::Int,nIncl::Int, FT)

    # First, get 8 scalars:
    n1 = ntuple(i->FT(0), 8)
    # get 3 Arrays (nlayer+1):
    n11 = ntuple(i->Array{Float32}(undef, nLayer+1), 3)
    # get 4 Arrays (nlayer)
    n12 = ntuple(i->Array{Float32}(undef, nLayer), 4)
    # get wl dimension variables:
    n2 = ntuple(i->Array{Float32}(undef, nWL), 6)
    # get wavelength and layers
    n3 = ntuple(i->Array{Float32}(undef,nLayer+1, nWL), 2)
    # get full angles (3D)
    n4 = ntuple(i->Array{Float32}(undef, nIncl, nAzi, nLayer), 4)
    return struct_canopyRadiation{FT}(n1...,n11...,n12..., n2..., n3..., n4...)
end

function loadOpti(swl::Vector; file=file_Opti)
    # Read in all optical data:
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

    nr = typeof(swl)(undef,length(swl)-1)
    Km = typeof(swl)(undef,length(swl)-1)
    Kab = typeof(swl)(undef,length(swl)-1)
    Kant = typeof(swl)(undef,length(swl)-1)
    Kcar = typeof(swl)(undef,length(swl)-1)
    Kw = typeof(swl)(undef,length(swl)-1)
    KBrown = typeof(swl)(undef,length(swl)-1)
    phi = typeof(swl)(undef,length(swl)-1)
    KcaV = typeof(swl)(undef,length(swl)-1)
    KcaZ = typeof(swl)(undef,length(swl)-1)
    lambda = typeof(swl)(undef,length(swl)-1)
    kChlrel = typeof(swl)(undef,length(swl)-1)
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


function loadSun(swl::Vector; file=file_Sun)
    # Read in all optical data:
    suni = matread(file_Sun)["sun"]
    wl   =  suni["wl"]
    Edir =  suni["Edirect"]
    Ediff =  suni["Ediffuse"]

    wl_ = typeof(swl)(undef,length(swl)-1)
    Edir_ = typeof(swl)(undef,length(swl)-1)
    Ediff_ = typeof(swl)(undef,length(swl)-1)
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
