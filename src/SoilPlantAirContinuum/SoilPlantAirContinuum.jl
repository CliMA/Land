module SoilPlantAirContinuum

using ConstrainedRootSolvers: ReduceStepMethodND, SolutionToleranceND, find_peak
using DocStringExtensions: TYPEDFIELDS
using EmeraldConstants: GAS_R, GRAVITY, K_STEFAN, K_VON_KARMAN, M_H₂O, P_ATM, RT₂₅, T₀, T₂₅, YEAR_D, ρ_H₂O
using PkgUtility: numerical∫, tinfo
using WaterPhysics: latent_heat_vapor, relative_diffusive_coefficient, saturation_vapor_pressure

using ..CanopyLayers: LAND_2017, LAND_2021
using ..CanopyLayers: Canopy4RT, CanopyOpticals, CanopyRads, IncomingRadiation, LeafBios, RTCache, RTDimensions, SIF_fluxes!, SoilOpticals, SolarAngles, WaveLengths,
      canopy_fluxes!, canopy_geometry!, canopy_matrices!, fluspect!, short_wave!, thermal_fluxes!
using ..Photosynthesis: AbstractPhotoModelParaSet, AirLayer, C3CLM, C3ParaSet, C4ParaSet, GCO₂Mode, Leaf, leaf_photosynthesis!, leaf_rd!, leaf_temperature_dependence!
using ..PlantHydraulics: GrassLikeOrganism, PalmLikeOrganism, SteadyStateMode, TreeLikeOrganism, TreeSimple, create_grass, critical_flow, end_pressure, flow_profile!,
      pressure_profile!, roots_flow!, soil_p_25_swc, soil_swc, temperature_effects!
using ..StomataModels: BetaGLinearPsoil, CanopyLayer, ESMBallBerry, ESMMedlyn, EmpiricalStomatalModel, GswDrive, OptimizationStomatalModel, gas_exchange!, gsw_control!, prognostic_gsw!,
      stomatal_conductance, update_leaf_TP!, β_factor


# define constants
# TODO add unit conversion in PkgUtility
KG_2_MOL(FT)     = 1 / M_H₂O(FT);
KG_H_2_MOL_S(FT) = KG_2_MOL(FT) / 3600;


include("types/spacmono.jl");

include("layers/beta.jl");
include("layers/layer_fluxes.jl");
include("layers/initializert.jl");
include("layers/test_diurnal.jl");
include("layers/test_soil.jl");
include("layers/windspeed.jl");

include("planet/atmpressure.jl");
include("planet/solarangle.jl");

include("quantities/measures.jl");
include("quantities/photosynthesis.jl");

include("wind.jl");


end # module
