module Plant

export Tree,
       initialize_rt_module,
       update_canopy_from_rt_module!,
       update_tree_e_crit!,
       update_tree_with_time!

# load modules and packages
using CLIMAParameters
using DocStringExtensions
using Parameters

# load sub-module CanopyRT
using ..CanopyRT
@unpack computeCanopyGeomProps!,
        computeCanopyMatrices!,
        computeThermalFluxes!,
        create_canopyOpt,
        deriveCanopyFluxes!,
        fluspect!,
        leafbio,
        RTM_SW!,
        struct_canopy,
        struct_canopyOptProps,
        struct_canopyRadiation = CanopyRT

using ..PhotosynthesisModels
@unpack AbstractPhotoModelParaSet,
        C3VcVpJBernacchi,
        C4VcVpJCLM,
        get_an_ag_r_pi_from_gsc,
        get_an_ag_r_pi_from_gsc_list = PhotosynthesisModels

# include the constants
include("constants.jl")

# include the tree types
include("types/tree_canopy.jl")
include("types/tree_root.jl"  )
include("types/tree_stem.jl"  )
include("types/tree.jl"       )

# include the dynamic functions
include("dynamic/update_tree_with_time.jl")

# include the hydraulic functions
include("hydraulics/get_p_base_q_list_from_q.jl")
include("hydraulics/get_q_layer_from_p_base.jl" )
include("hydraulics/get_struct_p_end_from_q.jl" )
include("hydraulics/update_struct_from_q.jl"    )
include("hydraulics/update_tree_e_crit.jl"      )

# include the interface functions
include("interface/initialize_rt_module.jl"        )
include("interface/update_canopy_from_rt_module.jl")

# include the stomatal optimization functions
include("stomata/get_marginal_gain.jl"        )
include("stomata/get_marginal_penalty_wang.jl")

# include the water property functions, all temperature in K by default
include("water/get_relative_surface_tension.jl")
include("water/get_relative_viscosity.jl"      )
include("water/get_saturated_vapor_pressure.jl")
include("water/get_specific_latent_heat.jl"    )

# include the testing functions and tools
include("tool/get_a_par_curve.jl")
include("tool/get_a_pi_curve.jl" )

end
