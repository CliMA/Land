module SoilPlantAirContinuum

using ..CanopyLayers
using CLIMAParameters
using CLIMAParameters.Planet
using CLIMAParameters.SubgridScale
using ConstrainedRootSolvers
using DataFrames
using DocStringExtensions
using Parameters
using ..Photosynthesis
using ..PlantHydraulics
using Statistics
using ..StomataModels
using WaterPhysics




# define constants
struct EarthParameterSet <: AbstractEarthParameterSet end
const EARTH       = EarthParameterSet();
GAS_R(FT)         = FT( gas_constant() );
GRAVITY(FT)       = FT( grav(EARTH) );
K_0(FT)           = FT( T_freeze(EARTH) );
K_25(FT)          = K_0(FT) + 25;
K_STEFAN(FT)      = FT( Stefan() );
K_VON_KARMAN(FT)  = FT( von_karman_const(EARTH) );
MOLMASS_WATER(FT) = FT( molmass_water(EARTH) );
P_ATM(FT)         = FT( MSLP(EARTH) );
RK_25(FT)         = GAS_R(FT) * K_25(FT);
YEAR_D(FT)        = FT( 365.2422222 );
ρ_H₂O(FT)         = FT( ρ_cloud_liq(EARTH) );
KG_2_MOL(FT)      = 1 / MOLMASS_WATER(FT);
KG_H_2_MOL_S(FT)  = KG_2_MOL(FT) / 3600;




# export public types
export SPACContainer1L,
       SPACContainer2L,
       SPACMono,
       SPACSimple




# export public functions
export annual_profit,
       annual_simulation!,
       atmospheric_pressure,
       atmospheric_pressure_ratio,
       big_leaf_partition!,
       create_dataframe,
       gain_risk_map,
       initialize_spac_canopy!,
       leaf_allocation!,
       leaf_gas_exchange!,
       leaf_gas_exchange_nonopt!,
       leaf_temperature,
       leaf_temperature_shaded,
       leaf_temperature_sunlit,
       optimize_flows!,
       optimize_hs!,
       optimize_leaf!,
       ppm_to_Pa,
       test_soil_from_psoil,
       test_soil_from_swc,
       vary_spac!,
       zenith_angle




include("types/container.jl" )
include("types/spacmono.jl"  )
include("types/spacsimple.jl")

include("bigleaf/gainriskmap.jl"   )
include("bigleaf/gasexchange.jl"   )
include("bigleaf/leafallocation.jl")
include("bigleaf/optimizeflow.jl"  )
include("bigleaf/optimizehs.jl"    )
include("bigleaf/optimizeleaf.jl"  )
include("bigleaf/partition.jl"     )
include("bigleaf/temperature.jl"   )
include("bigleaf/varytrait.jl"     )

include("layers/initializert.jl")
include("layers/test_diurnal.jl")
include("layers/test_soil.jl"   )
include("layers/windspeed.jl"   )

include("planet/atmpressure.jl")
include("planet/solarangle.jl" )

include("simulation/annualprofit.jl"    )
include("simulation/annualsimulation.jl")
include("simulation/createdataframe.jl" )

include("soil/moisture.jl")




end # module
