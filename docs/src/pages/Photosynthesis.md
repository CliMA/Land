# Photosynthesis module for the LAND model

## Temperature Dependency
### Types
```@docs
Land.Photosynthesis.AbstractTDParameterSet
Land.Photosynthesis.ArrheniusTD
Land.Photosynthesis.ArrheniusPeakTD
```

### Parameter sets
```@docs
Land.Photosynthesis.JmaxTDBernacchi
Land.Photosynthesis.JmaxTDCLM
Land.Photosynthesis.JmaxTDLeuning
Land.Photosynthesis.KcTDBernacchi
Land.Photosynthesis.KcTDCLM
Land.Photosynthesis.KoTDBernacchi
Land.Photosynthesis.KoTDCLM
Land.Photosynthesis.KpepTDBoyd
Land.Photosynthesis.KpepTDCLM
Land.Photosynthesis.RespirationTDBernacchi
Land.Photosynthesis.RespirationTDCLM
Land.Photosynthesis.VcmaxTDBernacchi
Land.Photosynthesis.VcmaxTDCLM
Land.Photosynthesis.VcmaxTDLeuning
Land.Photosynthesis.VomaxTDBernacchi
Land.Photosynthesis.VpmaxTDBoyd
Land.Photosynthesis.ΓStarTDBernacchi
Land.Photosynthesis.ΓStarTDCLM
```

## Photosynthesis Colimitation
```@docs
Land.Photosynthesis.AbstractColimitation
Land.Photosynthesis.MinColimit
Land.Photosynthesis.CurvedColimit
```

## Fluorescence Model
### Types
```@docs
Land.Photosynthesis.AbstractFluoModelParaSet
Land.Photosynthesis.FluoParaSet
```

### Parameter sets
```@docs
Land.Photosynthesis.FluorescenceFlexas
```

## Respiration Quantification
```@docs
Land.Photosynthesis.VtoRCollatz
Land.Photosynthesis.VtoRDefault
```

## Stomatal Model Scheme
```@docs
Land.Photosynthesis.AbstractStomatalModel
Land.Photosynthesis.EmpiricalStomatalModel
Land.Photosynthesis.OptimizationStomatalModel
```

### Empirical Stomatal Model
```@docs
Land.Photosynthesis.ESMBallBerry
Land.Photosynthesis.ESMGentine
Land.Photosynthesis.ESMLeuning
Land.Photosynthesis.ESMMedlyn
```

### Optimization Stomtal Model
```@docs
Land.Photosynthesis.OSMEller
Land.Photosynthesis.OSMSperry
Land.Photosynthesis.OSMWang
Land.Photosynthesis.OSMWAP
Land.Photosynthesis.OSMWAPMod
```

## Parameter Sets Collection
### Types
```@docs
Land.Photosynthesis.AbstractPhotoModelParaSet
Land.Photosynthesis.C3ParaSet
Land.Photosynthesis.C4ParaSet
```

### Parameter sets
```@docs
Land.Photosynthesis.C3Bernacchi
Land.Photosynthesis.C3CLM
Land.Photosynthesis.C4CLM
```

## Leaf and Environment
```@docs
Land.Photosynthesis.Leaf
Land.Photosynthesis.AirLayer
```

## Photosynthesis Model
### Temperature corrections
```@docs
Land.Photosynthesis.arrhenius_correction
Land.Photosynthesis.photo_TD_from_set
Land.Photosynthesis.photo_TD_from_val
Land.Photosynthesis.leaf_jmax!
Land.Photosynthesis.leaf_kc!
Land.Photosynthesis.leaf_km!
Land.Photosynthesis.leaf_ko!
Land.Photosynthesis.leaf_kpep!
Land.Photosynthesis.leaf_rd!
Land.Photosynthesis.leaf_vcmax!
Land.Photosynthesis.leaf_vpmax!
Land.Photosynthesis.leaf_Γstar!
```

### Photosynthetic rates (individual)
```@docs
Land.Photosynthesis.rubisco_limited_rate!
Land.Photosynthesis.leaf_ETR!
Land.Photosynthesis.light_limited_rate!
Land.Photosynthesis.product_limited_rate!
Land.Photosynthesis.photosynthesis_colimit
Land.Photosynthesis.leaf_fluorescence!
```

### Photosynthetic rates (collection)
```@docs
Land.Photosynthesis.photo_temperature_dependence!
Land.Photosynthesis.photo_radiation_dependence!
Land.Photosynthesis.photo_CO₂_dependence!
Land.Photosynthesis.leaf_photo_from_pi!
```

### Photosynthesis from Stomatal Conductance
```@docs
Land.Photosynthesis.glc_diff!
Land.Photosynthesis.leaf_photo_from_glc!
```

## Stomatal Response to Environment
```@docs
Land.Photosynthesis.empirical_gsw_from_model
Land.Photosynthesis.envir_diff!
Land.Photosynthesis.leaf_photo_from_envir!
Land.Photosynthesis.leaf_gsw_control!
Land.Photosynthesis.leaf_heat_flux!
```
