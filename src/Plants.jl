module Plants

using CLIMAParameters
using DataFrames
using DocStringExtensions
using Parameters
using Photosynthesis
using PlantHydraulics
using StomataModels
using WaterPhysics

Planet = CLIMAParameters.Planet




include("constants.jl")

include("types/container.jl")
include("types/spac.jl"     )

include("planet/atmpressure.jl")
include("planet/solarangle.jl" )

include("canopy/bigleaf.jl"    )
include("canopy/gasexchange.jl")
include("canopy/temperature.jl")




include("yujie_111/simulate_growing_season.jl")
include("./yujie_111/yujie_111_evaluate_model.jl"        )
include("./yujie_111/yujie_111_gain_risk_matrix.jl"      )
include("./yujie_111/yujie_111_get_annual_cica.jl"       )
include("./yujie_111/yujie_111_get_annual_profit.jl"     )
include("./yujie_111/yujie_111_get_optimal_bcklv.jl"     )
include("./yujie_111/yujie_111_get_optimal_dade.jl"      )
include("./yujie_111/yujie_111_get_optimal_f_knowt.jl"   )
include("./yujie_111/yujie_111_get_optimal_fs.jl"        )
include("./yujie_111/yujie_111_get_optimal_fs_map.jl"    )
include("./yujie_111/yujie_111_get_optimal_investment.jl")
include("./yujie_111/yujie_111_get_pacg_knowt.jl"        )
include("./yujie_111/yujie_111_get_pacgt.jl"             )
include("./yujie_111/yujie_111_get_pacgts.jl"            )
include("./yujie_111/yujie_111_set_soil_type.jl"         )
include("./yujie_111/yujie_111_update_leaf.jl"           )
include("./yujie_111/yujie_111_update_leaf_area.jl"      )
include("./yujie_111/yujie_111_update_soil.jl"           )
include("./yujie_111/yujie_111_update_soil_p.jl"         )
include("./yujie_111/yujie_111_update_soil_swc.jl"       )


end # module
