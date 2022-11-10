module LeafOptics

using SpecialFunctions: expint
using UnPack: @unpack

using EmeraldConstants: AVOGADRO, H_PLANCK, LIGHT_SPEED, M_H₂O, ρ_H₂O
using ClimaCache: HyperspectralAbsorption, HyperspectralRadiation, HyperspectralLeafBiophysics, WaveLengthSet
using ClimaCache: MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC


include("fluorescence.jl" )
include("photon.jl"       )
include("radiation.jl"    )
include("spectra.jl"      )
include("transmittance.jl")


end # module
