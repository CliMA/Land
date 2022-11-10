module StomataModels

using UnPack: @unpack

using ..EmeraldConstants: CP_D_MOL, K_STEFAN, M_H₂O
using ..WaterPhysics: latent_heat_vapor, relative_diffusive_coefficient, relative_surface_tension, relative_viscosity, saturation_vapor_pressure
using ..ClimaCache: AbstractSoilVC, AbstractXylemVC
using ..ClimaCache: AndereggSM, BallBerrySM, BetaParameterG1, BetaParameterVcmax, EllerSM, GentineSM, LeuningSM, MedlynSM, SperrySM, WangSM, Wang2SM
using ..ClimaCache: AirLayer, BroadbandLeafBiophysics, C4VJPModel, GCO₂Mode, HyperspectralLeafBiophysics, Leaf, LeafHydraulics, Leaves1D, Leaves2D
using ..ClimaCache: MonoElementSPAC, MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC
using ..Photosynthesis: leaf_photosynthesis!, ∂R∂T
using ..SoilHydraulics: relative_hydraulic_conductance
using ..PlantHydraulics: relative_hydraulic_conductance, ∂E∂P


# include files
include("../../packages/StomataModels.jl/src/conductance.jl")
include("../../packages/StomataModels.jl/src/empirical.jl"  )
include("../../packages/StomataModels.jl/src/limits.jl"     )
include("../../packages/StomataModels.jl/src/optimality.jl" )


end # module
