# PHOTOSYNTHESIS MODELS MODULE

## Temperature Dependency

### Jmax
```@docs
Land.PhotosynthesisModels.get_jmax
Land.PhotosynthesisModels.get_j
Land.PhotosynthesisModels.JmaxTDBernacchi
Land.PhotosynthesisModels.JmaxTDCLM
Land.PhotosynthesisModels.JmaxTDLeuning
```

### Kc
```@docs
Land.PhotosynthesisModels.get_kc
Land.PhotosynthesisModels.KcTDBernacchi
Land.PhotosynthesisModels.KcTDCLM
```

### Ko
```@docs
Land.PhotosynthesisModels.get_ko
Land.PhotosynthesisModels.KoTDBernacchi
Land.PhotosynthesisModels.KoTDCLM
```

### Kpep
```@docs
Land.PhotosynthesisModels.get_kpep
Land.PhotosynthesisModels.KpepTDBoyd
Land.PhotosynthesisModels.KpepTDCLM
```

### Respiration
```@docs
Land.PhotosynthesisModels.get_r
Land.PhotosynthesisModels.RespirationTDBernacchi
Land.PhotosynthesisModels.RespirationTDCLM
```

### Vmax
```@docs
Land.PhotosynthesisModels.get_vmax
Land.PhotosynthesisModels.VcmaxTDBernacchi
Land.PhotosynthesisModels.VcmaxTDCLM
Land.PhotosynthesisModels.VcmaxTDLeuning
Land.PhotosynthesisModels.VomaxTDBernacchi
Land.PhotosynthesisModels.VpmaxTDBoyd
```

### Γ*
```@docs
Land.PhotosynthesisModels.get_Γ_star
Land.PhotosynthesisModels.ΓStarTDBernacchi
Land.PhotosynthesisModels.ΓStarTDCLM
```

## Arrhenius Correction
```@docs
Land.PhotosynthesisModels.arrhenius_correction
```

## Parameter Sets

### Respiration~Vcmax
```@docs
Land.PhotosynthesisModels.VtoRCollatz
Land.PhotosynthesisModels.VtoRDefault
```

### Preset Photosynthesis Types
```@docs
Land.PhotosynthesisModels.C3VcJBernacchi
Land.PhotosynthesisModels.C3VcVpJBernacchi
Land.PhotosynthesisModels.C4VcVpJBoyd
Land.PhotosynthesisModels.C4VcVpJCLM
```

## Procedure

### Calculate Photosynthesis from CO2
```@docs
Land.PhotosynthesisModels.get_an_ag_r_from_pi
```

### Calculate Photosynthesis from gsc
```@docs
Land.PhotosynthesisModels.get_an_ag_r_pi_from_gsc
Land.PhotosynthesisModels.get_an_ag_r_pi_from_gsc_list
```
