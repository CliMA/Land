module StomataModels

using CLIMAParameters
using ConstrainedRootSolvers
using DocStringExtensions
using Parameters
using Photosynthesis
using PlantHydraulics
using WaterPhysics

Planet = CLIMAParameters.Planet




# Define constants
struct EarthParameterSet <: AbstractEarthParameterSet end
const EARTH = EarthParameterSet()
const K_25          = Planet.T_freeze(EARTH) + 25
const MOLMASS_WATER = Planet.molmass_water(EARTH)




include("types.jl"    )
include("refresh.jl"  )
include("empirical.jl")
include("solution.jl" )
include("stomata.jl"  )




end # module
