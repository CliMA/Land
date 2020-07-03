# Photosynthesis

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
Land.Photosynthesis.rubisco_limited_rate_glc!
Land.Photosynthesis.leaf_ETR!
Land.Photosynthesis.light_limited_rate!
Land.Photosynthesis.light_limited_rate_glc!
Land.Photosynthesis.product_limited_rate!
Land.Photosynthesis.product_limited_rate_glc!
Land.Photosynthesis.leaf_fluorescence!
```

### Photosynthetic rates
```@docs
Land.Photosynthesis.leaf_temperature_dependence!
Land.Photosynthesis.leaf_photo_from_pi!
Land.Photosynthesis.leaf_photo_from_glc!
```
