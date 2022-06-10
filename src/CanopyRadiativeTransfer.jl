module CanopyRadiativeTransfer

using ClimaCache: HyperspectralMLCanopy, HyperspectralRadiation, Leaf, Soil, SunSensorGeometry, VerhoefLIDF
using DocStringExtensions: METHODLIST
using LeafOptics: photon
using LinearAlgebra: mul!
using QuadGK: quadgk
using Statistics: mean
using UnPack: @unpack


include("clumping.jl"         )
include("coefficients.jl"     )
include("geometry.jl"         )
include("inclination_angle.jl")
include("radiation.jl"        )


end # module
