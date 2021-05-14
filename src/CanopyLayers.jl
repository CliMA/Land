module CanopyLayers

using LazyArtifacts

using DocStringExtensions: TYPEDFIELDS
using LinearAlgebra: mul!
using MAT: matread
using PkgUtility: AVOGADRO, H_PLANCK, K_STEFAN, LIGHT_SPEED, numerical∫
using Polynomials: Polynomial
using QuadGK: quadgk
using Statistics: mean
using UnPack: @unpack




# define the constants
const OPTI_2021 = artifact"land_model_spectrum" *
                  "/Optipar2021_ProspectPRO_CX.mat";
const OPTI_2017 = artifact"land_model_spectrum" * "/Optipar2017_ProspectD.mat";
const FILE_SUN  = artifact"land_model_spectrum" * "/sun.mat";




# export public types
export Canopy4RT, CanopyOpticals, CanopyRads, IncomingRadiation, LeafBios,
       LeafOpticals, RTCache, RTDimensions, SoilOpticals, SolarAngles,
       WaveLengths

# export public functions
export big_leaf_partition, canopy_fluxes!, canopy_geometry!, canopy_matrices!,
       create_canopy_opticals, create_canopy_rads, create_canopy_rt,
       create_incoming_radiation, create_leaf_bios, create_leaf_opticals,
       create_rt_cache, create_rt_dims, create_soil_opticals,
       create_wave_length, diffusive_S, fluspect!, initialize_rt_module,
       short_wave!, SIF_fluxes!, thermal_fluxes!

# Vegetation indices
export BLUE, EVI, EVI2, LSWI, NDVI, NIR, NIRv, NIRvES, RED, REF_WL, SIF_740,
       SIF_757, SIF_771, SIF_WL, SWIR




include("utils/calctav.jl"    )
include("utils/dladgen.jl"    )
include("utils/e2phot.jl"     )
include("utils/expint.jl"     )
include("utils/psofunction.jl")
include("utils/volscatt.jl"   )

include("types/canopy4rt.jl"     )
include("types/canopyopticals.jl")
include("types/canopyrads.jl"    )
include("types/caches.jl"        )
include("types/incomingrad.jl"   )
include("types/leafbios.jl"      )
include("types/leafopticals.jl"  )
include("types/rtdims.jl"        )
include("types/soilopticals.jl"  )
include("types/solarangles.jl"   )
include("types/wavelength.jl"    )

include("initialize/all.jl"           )
include("initialize/canopy4rt.jl"     )
include("initialize/canopyopticals.jl")
include("initialize/canopyrads.jl"    )
include("initialize/caches.jl"        )
include("initialize/incomingrad.jl"   )
include("initialize/leafbios.jl"      )
include("initialize/leafopticals.jl"  )
include("initialize/rtdims.jl"        )
include("initialize/soilopticals.jl"  )
include("initialize/wavelength.jl"    )

include("bigleaf/bigleaf.jl")

include("layers/canopyfluxes.jl"  )
include("layers/canopygeometry.jl")
include("layers/canopymatrices.jl")
include("layers/diffusives.jl"    )
include("layers/fluspect.jl"      )
include("layers/indicies.jl"      )
include("layers/shortwave.jl"     )
include("layers/siffluxes.jl"     )
include("layers/thermalfluxes.jl" )




end # module
