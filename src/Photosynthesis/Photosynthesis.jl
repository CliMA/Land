module Photosynthesis

using DocStringExtensions
using Parameters
using RootSolvers

using ..LandParameters
using ..MathTools
using ..RootSolversExtension
using ..WaterPhysics
using ..Hydraulics

# TODO Add β function options
# TODO update p_O₂ from p_H₂O
# TODO Add control for envir.p_H₂O > leaf.p_sat
# TODO Add a judgment of g_sw > g_max
# TODO A=Float32 does not work well

# define constants here
const CP_D             = LandParameters.CP_D
const GAS_R            = LandParameters.GAS_R
const K_25             = LandParameters.K_25
const MOLMASS_DRYAIR   = LandParameters.MOLMASS_DRYAIR
const MOLMASS_WATER    = LandParameters.MOLMASS_WATER

# export public types
export AirLayer,
       C3CLM,
       C4CLM,
       Leaf

# export public functions
export leaf_photo_from_envir!,
       leaf_photo_from_glc!,
       leaf_photo_from_pi!,
       leaf_temperature_dependence!




include("types.jl"     )
include("parasets.jl"  )
include("photomodel.jl")
include("stomata.jl"   )
include("thermo.jl"    )




end