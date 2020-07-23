module Plants

using CanopyRadiation
using CLIMAParameters
using ConstrainedRootSolvers
using DataFrames
using DocStringExtensions
using Parameters
using Photosynthesis
using PlantHydraulics
using Revise
using Statistics
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




# export public types
export SPACContainer1L,
       SPACContainer2L,
       SPACSimple




# export public functions
export annual_profit,
       annual_simulation!,
       atmospheric_pressure,
       atmospheric_pressure_ratio,
       big_leaf_partition!,
       gain_risk_map,
       leaf_allocation!,
       leaf_gas_exchange!,
       leaf_gas_exchange_nonopt!,
       leaf_temperature,
       leaf_temperature_shaded,
       leaf_temperature_sunlit,
       optimize_flows!,
       optimize_leaf!,
       ppm_to_Pa,
       zenith_angle




include("types/container.jl" )
include("types/spacmono.jl"  )
include("types/spacsimple.jl")

include("bigleaf/gainriskmap.jl" )
include("bigleaf/gasexchange.jl" )
include("bigleaf/optimizeflow.jl")
include("bigleaf/partition.jl"   )
include("bigleaf/temperature.jl" )

include("investment/leafallocation.jl")
include("investment/optimizeleaf.jl"  )

include("layers/test_diurnal.jl")

include("planet/atmpressure.jl")
include("planet/solarangle.jl" )

include("simulation/annualprofit.jl"    )
include("simulation/annualsimulation.jl")

include("soil/moisture.jl")




end # module
