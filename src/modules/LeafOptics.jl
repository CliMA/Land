module LeafOptics

using SpecialFunctions: expint
using UnPack: @unpack

using ..EmeraldConstants: AVOGADRO, H_PLANCK, LIGHT_SPEED, M_H₂O, ρ_H₂O
using ..EmeraldNamespace: HyperspectralAbsorption, HyperspectralRadiation, HyperspectralLeafBiophysics, WaveLengthSet
using ..EmeraldNamespace: MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC


include("../../packages/LeafOptics.jl/src/fluorescence.jl" )
include("../../packages/LeafOptics.jl/src/photon.jl"       )
include("../../packages/LeafOptics.jl/src/radiation.jl"    )
include("../../packages/LeafOptics.jl/src/spectra.jl"      )
include("../../packages/LeafOptics.jl/src/transmittance.jl")


end # module
