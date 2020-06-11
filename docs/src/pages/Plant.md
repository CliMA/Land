# Plant module in the LAND model

## Plant structure
```@docs
Land.Plant.RootLayer
Land.Plant.create_root_list
Land.Plant.Stem
Land.Plant.create_branch_list
Land.Plant.Leaf
Land.Plant.CanopyLayer
Land.Plant.Tree
```

## Stomatal models
```@docs
Land.Plant.AbstractStomatalModel
Land.Plant.AbstractEmpiricalStomatalModel
Land.Plant.ESMBallBerry
Land.Plant.ESMGentine
Land.Plant.ESMLeuning
Land.Plant.ESMMedlyn
Land.Plant.AbstractOptimizationStomatalModel
Land.Plant.OSMEller
Land.Plant.OSMSperry
Land.Plant.OSMWang
Land.Plant.OSMWAP
Land.Plant.OSMWAPMod
```

## Pant hydraulics
```@docs
Land.Plant.get_p_base_q_list_from_q
Land.Plant.get_q_layer_from_p_base
Land.Plant.get_struct_p_end_from_q
Land.Plant.update_struct_from_q!
Land.Plant.update_leaf_ak_max!
Land.Plant.update_tree_e_crit!
```

## Optimality stomatal model
```@docs
Land.Plant.get_marginal_gain
Land.Plant.get_marginal_penalty
```

## Empirical stomatal model
```@docs
Land.Plant.get_empirical_gsw_from_model
```

## Stomatal dynamics
```@docs
Land.Plant.update_tree_with_time!
```