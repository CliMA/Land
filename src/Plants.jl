module Plants

using CLIMAParameters
using DataFrames
using Parameters
using Photosynthesis
using PlantHydraulics
using Printf
using StomataModels
using WaterPhysics

Planet = CLIMAParameters.Planet

include("constants.jl")
include("Earth/Earth.jl")
include("Radiation/Radiation.jl")
include("Plant/Plant.jl")








include("./yujie_111/yujie_111.jl"                       )
include("./yujie_111/yujie_111_evaluate_model.jl"        )
include("./yujie_111/yujie_111_gain_risk_matrix.jl"      )
include("./yujie_111/yujie_111_get_annual_cica.jl"       )
include("./yujie_111/yujie_111_get_annual_gpp.jl"        )
include("./yujie_111/yujie_111_get_annual_profit.jl"     )
include("./yujie_111/yujie_111_get_e_crit.jl"            )
include("./yujie_111/yujie_111_get_k_tree.jl"            )
include("./yujie_111/yujie_111_get_leaf_partition.jl"    )
include("./yujie_111/yujie_111_get_optimal_bcklv.jl"     )
include("./yujie_111/yujie_111_get_optimal_dade.jl"      )
include("./yujie_111/yujie_111_get_optimal_f_knowt.jl"   )
include("./yujie_111/yujie_111_get_optimal_fs.jl"        )
include("./yujie_111/yujie_111_get_optimal_fs_map.jl"    )
include("./yujie_111/yujie_111_get_optimal_investment.jl")
include("./yujie_111/yujie_111_get_p.jl"                 )
include("./yujie_111/yujie_111_get_pacg_knowt.jl"        )
include("./yujie_111/yujie_111_get_pacgt.jl"             )
include("./yujie_111/yujie_111_get_pacgts.jl"            )
include("./yujie_111/yujie_111_get_pacgts_gpp.jl"        )
include("./yujie_111/yujie_111_get_ps.jl"                )
include("./yujie_111/yujie_111_get_t_leaf.jl"            )
include("./yujie_111/yujie_111_init_legacy.jl"           )
include("./yujie_111/yujie_111_set_soil_type.jl"         )
include("./yujie_111/yujie_111_update_legacy.jl"         )
include("./yujie_111/yujie_111_update_leaf.jl"           )
include("./yujie_111/yujie_111_update_leaf_area.jl"      )
include("./yujie_111/yujie_111_update_soil.jl"           )
include("./yujie_111/yujie_111_update_soil_p.jl"         )
include("./yujie_111/yujie_111_update_soil_swc.jl"       )


include("./simulate_growing_season.jl")


end # module
