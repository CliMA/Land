# Photosynthesis

## Temperature Dependency
### Types
```@docs
Photosynthesis.AbstractTDParameterSet
Photosynthesis.ArrheniusTD
Photosynthesis.ArrheniusPeakTD
```

### Parameter sets
```@docs
Photosynthesis.JmaxTDBernacchi
Photosynthesis.JmaxTDCLM
Photosynthesis.JmaxTDLeuning
Photosynthesis.KcTDBernacchi
Photosynthesis.KcTDCLM
Photosynthesis.KoTDBernacchi
Photosynthesis.KoTDCLM
Photosynthesis.KpepTDBoyd
Photosynthesis.KpepTDCLM
Photosynthesis.RespirationTDBernacchi
Photosynthesis.RespirationTDCLM
Photosynthesis.VcmaxTDBernacchi
Photosynthesis.VcmaxTDCLM
Photosynthesis.VcmaxTDLeuning
Photosynthesis.VomaxTDBernacchi
Photosynthesis.VpmaxTDBoyd
Photosynthesis.ΓStarTDBernacchi
Photosynthesis.ΓStarTDCLM
```

## Fluorescence Model
### Types
```@docs
Photosynthesis.AbstractFluoModelParaSet
Photosynthesis.FluoParaSet
```

### Parameter sets
```@docs
Photosynthesis.FluorescenceFlexas
```

## Respiration Quantification
```@docs
Photosynthesis.VtoRCollatz
Photosynthesis.VtoRDefault
```

## Parameter Sets Collection
### Types
```@docs
Photosynthesis.AbstractPhotoModelParaSet
Photosynthesis.C3ParaSet
Photosynthesis.C4ParaSet
```

### Parameter sets
```@docs
Photosynthesis.C3Bernacchi
Photosynthesis.C3CLM
Photosynthesis.C4CLM
```

## Leaf and Environment
```@docs
Photosynthesis.Leaf
Photosynthesis.AirLayer
```

## Photosynthesis Model
### Temperature corrections
```@docs
Photosynthesis.arrhenius_correction
Photosynthesis.photo_TD_from_set
Photosynthesis.photo_TD_from_val
Photosynthesis.leaf_jmax!
Photosynthesis.leaf_kc!
Photosynthesis.leaf_km!
Photosynthesis.leaf_ko!
Photosynthesis.leaf_kpep!
Photosynthesis.leaf_rd!
Photosynthesis.leaf_vcmax!
Photosynthesis.leaf_vpmax!
Photosynthesis.leaf_Γstar!
```

### Photosynthetic rates (individual)
```@docs
Photosynthesis.rubisco_limited_rate!
Photosynthesis.rubisco_limited_rate_glc!
Photosynthesis.leaf_ETR!
Photosynthesis.light_limited_rate!
Photosynthesis.light_limited_rate_glc!
Photosynthesis.product_limited_rate!
Photosynthesis.product_limited_rate_glc!
Photosynthesis.leaf_fluorescence!
```

### Photosynthetic rates
```@docs
Photosynthesis.leaf_temperature_dependence!
Photosynthesis.leaf_photo_from_pi!
Photosynthesis.leaf_photo_from_glc!
```
