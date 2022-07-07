module StomataModels

using ClimaCache: AbstractSoilVC, AbstractXylemVC, AirLayer, BallBerrySM, BetaFunction, BetaParameterG1, BetaParameterVcmax, C4VJPModel, GentineSM, Leaf, LeafHydraulics, Leaves1D, Leaves2D,
      LeuningSM, MedlynSM
using DocStringExtensions: METHODLIST
using PlantHydraulics: relative_hydraulic_conductance
using SoilHydraulics: relative_hydraulic_conductance
using UnPack: @unpack
using WaterPhysics: relative_diffusive_coefficient, relative_surface_tension


# include files
include("beta.jl"       )
include("conductance.jl")
include("empirical.jl"  )
include("limits.jl"     )


end # module
