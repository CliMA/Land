module SoilPlantAirContinuum

using CanopyRadiativeTransfer: canopy_fluorescence!, canopy_radiation!, soil_albedo!
using ClimaCache: AirLayer, GCO₂Mode, MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC
using EmeraldConstants: CP_L, CP_L_MOL, ρ_H₂O
using LeafOptics: leaf_spectra!
using Photosynthesis: leaf_photosynthesis!
using PlantHydraulics: flow_out, plant_energy!, xylem_flow_profile!, xylem_pressure_profile!
using SoilHydraulics: soil_budget!
using Statistics: mean
using StomataModels: stomatal_conductance!
using UnPack: @unpack
using WaterPhysics: saturation_vapor_pressure


include("budget.jl"    )
include("initialize.jl")
include("model.jl"     )
include("quantity.jl"  )
include("update.jl"    )


end # module
