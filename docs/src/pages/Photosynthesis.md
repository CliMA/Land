# Photosynthesis module for the LAND model

## Types and parameter sets

### Temperature dependency
```@docs
Land.Photosynthesis.AbstractTDParameterSet
Land.Photosynthesis.AbstractArrheniusTDParameterSet
Land.Photosynthesis.AbstractArrheniusPeakTDParameterSet
Land.Photosynthesis.ArrheniusTD
Land.Photosynthesis.ArrheniusPeakTD
Land.Photosynthesis.KcTDBernacchi
Land.Photosynthesis.KoTDBernacchi
Land.Photosynthesis.RespirationTDBernacchi
Land.Photosynthesis.VcmaxTDBernacchi
Land.Photosynthesis.VomaxTDBernacchi
Land.Photosynthesis.ΓStarTDBernacchi
Land.Photosynthesis.KpepTDBoyd
Land.Photosynthesis.VpmaxTDBoyd
Land.Photosynthesis.JmaxTDLeuning
Land.Photosynthesis.VcmaxTDLeuning
Land.Photosynthesis.KcTDCLM
Land.Photosynthesis.KoTDCLM
Land.Photosynthesis.KpepTDCLM
Land.Photosynthesis.ΓStarTDCLM
Land.Photosynthesis.JmaxTDBernacchi
Land.Photosynthesis.JmaxTDCLM
Land.Photosynthesis.RespirationTDCLM
Land.Photosynthesis.VcmaxTDCLM
Land.Photosynthesis.VtoRCollatz
Land.Photosynthesis.VtoRDefault
```

### Photosynthesis model
```@docs
Land.Photosynthesis.AbstractPhotoModelParaSet
Land.Photosynthesis.C3ParaSet
Land.Photosynthesis.C4ParaSet
Land.Photosynthesis.C3Bernacchi
Land.Photosynthesis.C3CLM
Land.Photosynthesis.C4CLM
```

### Fluorescence
```@docs
Land.Photosynthesis.AbstractFluoModelParaSet
Land.Photosynthesis.FluoParaSet
Land.Photosynthesis.FluorescenceFlexas
```

## Photosynthesis model
```@docs
Land.Photosynthesis.arrhenius_correction
Land.Photosynthesis.get_jmax
Land.Photosynthesis.get_j
Land.Photosynthesis.get_kc
Land.Photosynthesis.get_ko
Land.Photosynthesis.get_kpep
Land.Photosynthesis.get_r
Land.Photosynthesis.get_vmax
Land.Photosynthesis.get_Γ_star
Land.Photosynthesis.an_ag_r_from_pi
Land.Photosynthesis.an_ag_r_pi_from_gsc
```
