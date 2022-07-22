module LeafOptics

using ClimaCache: HyperspectralAbsorption, HyperspectralRadiation, HyperspectralLeafBiophysics, MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC, WaveLengthSet
using ClimaCache: AVOGADRO, H_PLANCK, LIGHT_SPEED, M_H₂O, ρ_H₂O
using SpecialFunctions: expint
using UnPack: @unpack


# export public functions
export leaf_PAR, leaf_SIF, leaf_spectra!


include("fluorescence.jl" )
include("photon.jl"       )
include("radiation.jl"    )
include("spectra.jl"      )
include("transmittance.jl")


end # module
