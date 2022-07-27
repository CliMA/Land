module SoilPlantAirContinuum

using CanopyRadiativeTransfer: canopy_radiation!
using ClimaCache: AirLayer, GCO₂Mode, MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC
using ClimaCache: CP_L, CP_L_MOL, ρ_H₂O
using LeafOptics: leaf_spectra!
using Photosynthesis: leaf_photosynthesis!
using PlantHydraulics: plant_energy!, xylem_flow_profile!, xylem_pressure_profile!
using SoilHydraulics: soil_budget!
using StomataModels: stomatal_conductance!
using UnPack: @unpack
using WaterPhysics: saturation_vapor_pressure


include("budget.jl"    )
include("initialize.jl")
include("model.jl"     )
include("update.jl"    )


end # module
