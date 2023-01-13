module LeafOptics

using SpecialFunctions: expint

using ..EmeraldConstants: M_H₂O, ρ_H₂O
using ..EmeraldOptics: average_transmittance, energy, energy!, photon, photon!
using ..EmeraldNamespace: HyperspectralAbsorption, HyperspectralRadiation, HyperspectralLeafBiophysics, WaveLengthSet
using ..EmeraldNamespace: MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC


include("../../packages/LeafOptics.jl/src/fluorescence.jl")
include("../../packages/LeafOptics.jl/src/radiation.jl"   )
include("../../packages/LeafOptics.jl/src/spectra.jl"     )


end # module
