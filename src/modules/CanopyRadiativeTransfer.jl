module CanopyRadiativeTransfer

using ConstrainedRootSolvers: ReduceStepMethodND, SolutionToleranceND, find_peak
using LinearAlgebra: mul!, pinv
using QuadGK: quadgk
using Statistics: mean
using UnPack: @unpack

using ..EmeraldConstants: K_STEFAN
using ..ClimaCache: BroadbandRadiation, BroadbandSLCanopy, BroadbandSoilAlbedo
using ..ClimaCache: HyperspectralMLCanopy, HyperspectralRadiation, HyperspectralSoilAlbedo
using ..ClimaCache: Leaves1D, Leaves2D, Soil, SunSensorGeometry, VerhoefLIDF
using ..ClimaCache: MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC
using ..LeafOptics: energy!, photon, photon!


include("../../packages/CanopyRadiativeTransfer.jl/src/constants.jl"        )
include("../../packages/CanopyRadiativeTransfer.jl/src/clumping.jl"         )
include("../../packages/CanopyRadiativeTransfer.jl/src/coefficients.jl"     )
include("../../packages/CanopyRadiativeTransfer.jl/src/fluorescence.jl"     )
include("../../packages/CanopyRadiativeTransfer.jl/src/geometry.jl"         )
include("../../packages/CanopyRadiativeTransfer.jl/src/inclination_angle.jl")
include("../../packages/CanopyRadiativeTransfer.jl/src/radiation.jl"        )
include("../../packages/CanopyRadiativeTransfer.jl/src/remote_sensing.jl"   )
include("../../packages/CanopyRadiativeTransfer.jl/src/soil.jl"             )


end # module
