# StomtaModels API
```@meta
CurrentModule = StomataModels
```

## Stomatal model schemes
```@docs
AbstractStomatalModel
EmpiricalStomatalModel
ESMBallBerry
ESMGentine
ESMLeuning
ESMMedlyn
OptimizationStomatalModel
OSMEller
OSMSperry
OSMWang
OSMWAP
OSMWAPMod
```

## Leaves
```@docs
CanopyLayer
```

## Stomatal conductance
```@docs
empirical_gsw_from_model
envir_diff!
leaf_photo_from_envir!
```

## Photosynthesis
```@docs
update_leaf_TP!
update_leaf_AK!
update_leaf_from_glc!
update_leaf_from_gsw!
leaf_gsw_control!
```
