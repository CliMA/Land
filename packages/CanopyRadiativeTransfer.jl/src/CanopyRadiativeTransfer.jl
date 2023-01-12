module CanopyRadiativeTransfer

using ConstrainedRootSolvers: ReduceStepMethodND, SolutionToleranceND, find_peak
using LinearAlgebra: mul!, pinv
using QuadGK: quadgk
using Statistics: mean
using UnPack: @unpack

#using EmeraldConstants: K_STEFAN
#using EmeraldNamespace: BroadbandRadiation, BroadbandSLCanopy, BroadbandSoilAlbedo
#using EmeraldNamespace: HyperspectralMLCanopy, HyperspectralRadiation, HyperspectralSoilAlbedo
#using EmeraldNamespace: Leaves1D, Leaves2D, Soil, SunSensorGeometry, VerhoefLIDF
#using EmeraldNamespace: MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC
#using LeafOptics: energy!, photon, photon!


include("constants.jl"        )
include("clumping.jl"         )
include("coefficients.jl"     )
include("fluorescence.jl"     )
include("geometry.jl"         )
include("inclination_angle.jl")
include("radiation.jl"        )
include("remote_sensing.jl"   )
include("soil.jl"             )


end # module
