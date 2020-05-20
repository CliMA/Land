# Land model

# Example Documentation

## Canopy RT module

### Leaf Bio Structure:
```@docs
Land.CanopyRT.leafbio
```

### Canopy Structure
```@docs
Land.CanopyRT.struct_canopy
```

### Canopy Radiation Structure
```@docs
Land.CanopyRT.struct_canopyRadiation
```

### Canopy Optical Property Structure
```@docs
Land.CanopyRT.struct_canopyOptProps
```

## Photosynthesis module

### Leaf Params:
```@docs
Land.Photosynthesis.leaf_params
```
### Leaf Rate Limiting Steps:
```@docs
Land.Photosynthesis.max_carboxylation_rate!
Land.Photosynthesis.max_electron_transport_rate!
Land.Photosynthesis.michaelis_menten_constants!
Land.Photosynthesis.leaf_respiration!
```


## Plant module

### Plant structure
```@docs
Land.Plant.Tree
Land.Plant.Root
Land.Plant.RootLayer
Land.Plant.Stem
Land.Plant.Branch
Land.Plant.Canopy
Land.Plant.CanopyLayer
Land.Plant.Leaf
```

### Plant hydraulics
```@docs
Land.Plant.get_struct_p_end_from_q
Land.Plant.get_q_layer_from_p_base
Land.Plant.get_p_base_q_list_from_q
Land.Plant.update_struct_from_q!
Land.Plant.update_tree_e_crit!
```

### Stoma dynamics
```@docs
Land.Plant.get_marginal_gain
Land.Plant.get_marginal_penalty_wang
Land.Plant.update_tree_with_time!
```

### Interfaces
```@docs
Land.Plant.update_canopy_from_rt_module!
```
