# directly imclude this file will save the time to compile other Land modules

# load modules and packages
using CLIMAParameters
using DocStringExtensions
using Parameters

# include the constants
include("constants.jl")

# include the tree struct
include("struct/struct_tree_canopy.jl")
include("struct/struct_tree_root.jl"  )
include("struct/struct_tree_stem.jl"  )
include("struct/struct_tree.jl"       )

# include the dynamic functions
include("dynamic/update_tree_with_time.jl")

# include the hydraulic functions
include("hydraulics/get_p_base_q_list_from_q.jl")
include("hydraulics/get_q_layer_from_p_base.jl" )
include("hydraulics/get_struct_p_end_from_q.jl" )
include("hydraulics/update_struct_from_q.jl"    )
include("hydraulics/update_tree_e_crit.jl"      )

# include the interface functions
include("interface/update_canopy_from_rt_module.jl")

# include the photosynthesis functions
include("photosynthesis/get_leaf_an_ag_r_from_pi.jl"    )
include("photosynthesis/get_leaf_an_ag_r_pi_from_gsc.jl")
include("photosynthesis/get_leaf_j.jl"                  )
include("photosynthesis/get_leaf_jmax.jl"               )
include("photosynthesis/get_leaf_r_from_r25.jl"         )
include("photosynthesis/get_leaf_r_from_v25.jl"         )
include("photosynthesis/get_leaf_vcmax.jl"              )

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
