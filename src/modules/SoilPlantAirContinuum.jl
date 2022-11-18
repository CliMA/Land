module SoilPlantAirContinuum

using Statistics: mean
using UnPack: @unpack

using ..EmeraldConstants: CP_L, CP_L_MOL, T₀, ρ_H₂O
using ..WaterPhysics: saturation_vapor_pressure
using ..ClimaCache: AirLayer, GCO₂Mode, MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC
using ..LeafOptics: leaf_spectra!
using ..CanopyRadiativeTransfer: canopy_fluorescence!, canopy_radiation!, soil_albedo!
using ..Photosynthesis: leaf_photosynthesis!
using ..SoilHydraulics: soil_budget!
using ..PlantHydraulics: flow_out, plant_energy!, xylem_flow_profile!, xylem_pressure_profile!, β_factor
using ..StomataModels: stomatal_conductance!


include("../../packages/SoilPlantAirContinuum.jl/src/budget.jl"    )
include("../../packages/SoilPlantAirContinuum.jl/src/initialize.jl")
include("../../packages/SoilPlantAirContinuum.jl/src/model.jl"     )
include("../../packages/SoilPlantAirContinuum.jl/src/quantity.jl"  )
include("../../packages/SoilPlantAirContinuum.jl/src/update.jl"    )


end # module
