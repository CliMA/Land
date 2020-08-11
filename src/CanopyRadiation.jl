module CanopyRadiation

using CLIMAParameters
using DocStringExtensions
using LinearAlgebra
using MAT
using Parameters
using Polynomials
using QuadGK
using Revise
using Statistics




# define the constants
AVOGADRO(FT)    = FT( avogad() );
H_PLANCK(FT)    = FT( h_Planck() );
K_BOLTZMANN(FT) = FT( k_Boltzmann() );
LIGHT_SPEED(FT) = FT( light_speed() );
const file_Opti = joinpath(@__DIR__, "../data/Optipar2017_ProspectD.mat")
const file_Sun  = joinpath(@__DIR__, "../data/sun.mat")




# export public types
export Canopy4RT,
       CanopyOpticals,
       CanopyRads,
       IncomingRadiation,
       LeafBios,
       LeafOpticals,
       RTContainer,
       SoilOpticals,
       SolarAngles,
       WaveLengths




# export public functions
export big_leaf_partition,
       canopy_fluxes!,
       canopy_geometry!,
       canopy_matrices!,
       create_canopy_opticals,
       create_incoming_radiation,
       create_leaf_bios,
       create_leaf_opticals,
       create_rt_container,
       create_soil_opticals,
       create_wave_length,
       diffusive_S,
       fluspect!,
       initialize_rt_module,
       short_wave!,
       sif_fluxes!,
       thermal_fluxes!




include("utils/calctav.jl"    )
include("utils/dladgen.jl"    )
include("utils/e2phot.jl"     )
include("utils/expint.jl"     )
include("utils/integral.jl"   )
include("utils/psofunction.jl")
include("utils/volscatt.jl"   )

include("types/canopy4rt.jl"     )
include("types/canopyopticals.jl")
include("types/canopyrads.jl"    )
include("types/container.jl"     )
include("types/incomingrad.jl"   )
include("types/leafbios.jl"      )
include("types/leafopticals.jl"  )
include("types/soilopticals.jl"  )
include("types/solarangles.jl"   )
include("types/wavelength.jl"    )

include("initialize/canopyopticals.jl")
include("initialize/container.jl"     )
include("initialize/incomingrad.jl"   )
include("initialize/leafbios.jl"      )
include("initialize/leafopticals.jl"  )
include("initialize/soilopticals.jl"  )
include("initialize/wavelength.jl"    )

include("bigleaf/bigleaf.jl")

include("layers/canopyfluxes.jl"  )
include("layers/canopygeometry.jl")
include("layers/canopymatrices.jl")
include("layers/diffusives.jl"    )
include("layers/fluspect.jl"      )
include("layers/initialize.jl"    )
include("layers/shortwave.jl"     )
include("layers/siffluxes.jl"     )
include("layers/thermalfluxes.jl" )




end # module
