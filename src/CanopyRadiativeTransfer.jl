module CanopyRadiativeTransfer

using ClimaCache: BroadbandSoilAlbedo, HyperspectralMLCanopy, HyperspectralRadiation, HyperspectralSoilAlbedo, Leaf, Soil, SunSensorGeometry, VerhoefLIDF
using DocStringExtensions: METHODLIST
using LeafOptics: energy!, photon, photon!
using LinearAlgebra: mul!
using PkgUtility: K_STEFAN
using QuadGK: quadgk
using Statistics: mean
using UnPack: @unpack


const MODIS_BAND_1 = [ 620,  670];  # RED
const MODIS_BAND_2 = [ 841,  876];  # NIR
const MODIS_BAND_3 = [ 459,  479];  # BLUE
const MODIS_BAND_7 = [2105, 2155];  # SWIR

const OCO2_SIF_759 = [758.17, 759.20];  # SIF 757
const OCO2_SIF_770 = [769.62, 770.28];  # SIF 770
const OCO3_SIF_759 = [758.26, 759.28];  # SIF 757
const OCO3_SIF_770 = [769.67, 770.34];  # SIF 770

const TROPOMI_SIF_683 = [680, 685];     # SIF 683
const TROPOMI_SIF_747 = [735, 758];     # SIF 747
const TROPOMI_SIF_751 = [743, 758];     # SIF 751


include("clumping.jl"         )
include("coefficients.jl"     )
include("fluorescence.jl"     )
include("geometry.jl"         )
include("inclination_angle.jl")
include("radiation.jl"        )
include("remote_sensing.jl"   )
include("soil.jl"             )


end # module
