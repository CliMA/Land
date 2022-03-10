module SoilPlantAirContinuum

using CanopyLayers: Canopy4RT, CanopyOpticals, CanopyRads, IncomingRadiation, LeafBios, RTCache, RTDimensions, SIF_fluxes!, SoilOpticals, SolarAngles, WaveLengths, big_leaf_partition, canopy_fluxes!,
      canopy_geometry!, canopy_matrices!, create_canopy_opticals, create_canopy_rads, create_canopy_rt, create_incoming_radiation, create_leaf_bios, create_rt_cache, create_rt_dims,
      create_wave_length, fluspect!, short_wave!, thermal_fluxes!
using ClimaCache: AirLayer, GCO₂Mode, Leaf
using ConstrainedRootSolvers: ReduceStepMethodND, SolutionToleranceND, find_peak
using DataFrames: DataFrame
using DocStringExtensions: TYPEDFIELDS
using Photosynthesis: leaf_photosynthesis!, photosystem_temperature_dependence!
using PkgUtility: GAS_R, GRAVITY, K_STEFAN, K_VON_KARMAN, M_H₂O, P_ATM, RT_25, T_0, T_25, YEAR_D, tinfo, ρ_H₂O
using PlantHydraulics: AbstractPlantOrganism, GrassLikeOrganism, PalmLikeOrganism, SteadyStateMode, TreeLikeOrganism, TreeSimple, create_grass, critical_flow, end_pressure, flow_profile!,
      pressure_profile!, roots_flow!, soil_p_25_swc, soil_swc, temperature_effects!
using StomataModels: AbstractStomatalModel, CanopyLayer, ESMBallBerry, EmpiricalStomatalModel, GswDrive, gas_exchange!, gsw_control!, prognostic_gsw!, stomatal_conductance, update_leaf_TP!
using UnPack: @unpack
using WaterPhysics: latent_heat_vapor, relative_diffusive_coefficient, saturation_vapor_pressure




# define constants
# TODO add unit conversion in PkgUtility
KG_2_MOL(FT)     = 1 / M_H₂O(FT);
KG_H_2_MOL_S(FT) = KG_2_MOL(FT) / 3600;




# export public types
export SPACMono, SPACSimple

# export public functions
export annual_profit, annual_simulation!, atmospheric_pressure, atmospheric_pressure_ratio, big_leaf_partition!, create_dataframe, gain_risk_map, initialize_spac_canopy!, layer_fluxes!,
       leaf_allocation!, leaf_gas_exchange!, leaf_gas_exchange_nonopt!, leaf_temperature, leaf_temperature_shaded, leaf_temperature_sunlit, optimize_flows!, optimize_hs!, optimize_leaf!, ppm_to_Pa,
       test_soil_from_psoil!, test_soil_from_swc!, update_Cab!, update_Kmax!, update_LAI!, update_VJR!, update_VJRWW!, update_Weibull!, vary_spac!, zenith_angle




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

include("layers/layer_fluxes.jl")
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
