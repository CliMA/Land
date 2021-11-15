###############################################################################
#
# Canopy optical variables
#
###############################################################################
"""
    mutable struct CanopyOpticals{FT}

A struct for canopy optical properties

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct CanopyOpticals{FT}
    # local storage of dimension information
    "Number of azimuth angles"
    nAzi  ::Int = 36
    "Number of inclination agles"
    nIncl ::Int = 9
    "Number of canopy layers"
    nLayer::Int = 5
    "Number of wave lengths"
    nWL   ::Int = 10

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
    Ps ::Array{FT,1} = zeros(FT, nLayer+1)
    "Probability of directly viewing a leaf in viewing direction"
    Po ::Array{FT,1} = zeros(FT, nLayer+1)
    "Bi-directional probability of directly viewing a leaf (solar->canopy->viewing)"
    Pso::Array{FT,1} = zeros(FT, nLayer+1)

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
    "abs(fo)"
    absfo  ::Array{FT,2} = zeros(FT, (nIncl, nAzi))
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
    R_sd::Array{FT,2} = zeros(FT, (nWL, nLayer+1))
    "Effective layer reflectance (diffuse->diffuse)"
    R_dd::Array{FT,2} = zeros(FT, (nWL, nLayer+1))
    "Solar direct radiation per layer)"
    Es_::Array{FT,2}  = zeros(FT, (nWL, nLayer+1))
end
