module Plants

using CLIMAParameters
using DataFrames
using DocStringExtensions
using Optim
using Parameters
using Photosynthesis
using PlantHydraulics
using StomataModels
using WaterPhysics

Planet = CLIMAParameters.Planet




# define constants
struct EarthParameterSet <: AbstractEarthParameterSet end
const EARTH         = EarthParameterSet()
const GAS_R         = gas_constant()
const GRAVITY       = Planet.grav(EARTH)
const K_0           = Planet.T_freeze(EARTH)
const K_25          = K_0 + 25
const K_BOLTZMANN   = k_Boltzmann()
const MOLMASS_WATER = Planet.molmass_water(EARTH)
const P_ATM         = Planet.MSLP(EARTH)
const RK_25         = GAS_R * K_25
const YEAR_D        = 365.2422222
const ρ_H₂O         = Planet.ρ_cloud_liq(EARTH)

const KG_2_MOL      = 1 / MOLMASS_WATER
const KG_H_2_MOL_S  = KG_2_MOL / 3600




include("types/container.jl")
include("types/spac.jl"     )

include("canopy/bigleaf.jl"     )
include("canopy/gainriskmap.jl" )
include("canopy/gasexchange.jl" )
include("canopy/optimizeflow.jl")
include("canopy/temperature.jl" )
include("canopy/testing.jl"     )

include("investment/leafallocation.jl")
include("investment/optimizeleaf.jl"  )

include("planet/atmpressure.jl")
include("planet/solarangle.jl" )

include("simulation/annualprofit.jl"    )
include("simulation/annualsimulation.jl")

include("soil/moisture.jl")




end # module
