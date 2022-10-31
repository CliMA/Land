module StomataModels

using ClimaCache: AbstractSoilVC, AbstractXylemVC, AirLayer, AndereggSM, BallBerrySM, BetaParameterG1, BetaParameterVcmax, BroadbandLeafBiophysics, C4VJPModel, EllerSM, GCO₂Mode, GentineSM,
      HyperspectralLeafBiophysics, Leaf, LeafHydraulics, Leaves1D, Leaves2D, LeuningSM, MedlynSM, MonoElementSPAC, MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC, SperrySM, WangSM, Wang2SM
using EmeraldConstants: CP_D_MOL, K_STEFAN, M_H₂O
using Photosynthesis: leaf_photosynthesis!, ∂R∂T
using PlantHydraulics: relative_hydraulic_conductance, ∂E∂P
using SoilHydraulics: relative_hydraulic_conductance
using UnPack: @unpack
using WaterPhysics: latent_heat_vapor, relative_diffusive_coefficient, relative_surface_tension, relative_viscosity, saturation_vapor_pressure


# include files
include("conductance.jl")
include("empirical.jl"  )
include("limits.jl"     )
include("optimality.jl" )


end # module
