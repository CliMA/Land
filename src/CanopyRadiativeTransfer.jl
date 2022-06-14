module CanopyRadiativeTransfer

using ClimaCache: HyperspectralMLCanopy, HyperspectralRadiation, Leaf, Soil, SunSensorGeometry, VerhoefLIDF
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


include("clumping.jl"         )
include("coefficients.jl"     )
include("fluorescence.jl"     )
include("geometry.jl"         )
include("inclination_angle.jl")
include("radiation.jl"        )
include("remote_sensing.jl"   )


end # module
