module LeafOptics

using ClimaCache: HyperspectralAbsorption, HyperspectralRadiation, LeafBiophysics, WaveLengthSet
using DocStringExtensions: METHODLIST
using PkgUtility: H_PLANCK, LIGHT_SPEED, AVOGADRO, numerical∫
using SpecialFunctions: expint
using UnPack: @unpack


# export public types from ClimaCache
export HyperspectralAbsorption, HyperspectralRadiation, LeafBiophysics, WaveLengthSet

# export public functions
export leaf_PAR, leaf_SIF, leaf_spectra!


include("fluorescence.jl" )
include("photon.jl"       )
include("radiation.jl"    )
include("spectra.jl"      )
include("transmittance.jl")


end # module
