# CanopyRT module for the LAND model

## Canopy RT functions
These are the functions being used in the Canopy RT module to do all the fun stuff (need to add more text here later)
```@docs
Land.CanopyRT.compute_canopy_geometry!
Land.CanopyRT.compute_canopy_matrices!
Land.CanopyRT.simulate_short_wave!
Land.CanopyRT.derive_canopy_fluxes!
Land.CanopyRT.compute_diffusive_S
Land.CanopyRT.compute_thermal_fluxes!
Land.CanopyRT.initialize_rt_module
Land.CanopyRT.computeSIF_Fluxes!
Land.CanopyRT.fluspect!
Land.CanopyRT.calctav
```


## Canopy RT structures

### Canopy
```@docs
Land.CanopyRT.Canopy4RT
```

### Canopy optical parameters
```@docs
Land.CanopyRT.AbstractCanopyOpti
Land.CanopyRT.CanopyOptiArray
Land.CanopyRT.CanopyOptiMArray
Land.CanopyRT.create_canopy_optical
```

### Canopy radiation
```@docs
Land.CanopyRT.CanopyRadiation
Land.CanopyRT.AbstractIncomingRadiation
Land.CanopyRT.IncomingRadiationArray
Land.CanopyRT.IncomingRadiationMArray
Land.CanopyRT.create_incoming_radiation
```

### Leaf biological parameters
```@docs
Land.CanopyRT.AbstractLeafBio
Land.CanopyRT.LeafBioArray
Land.CanopyRT.LeafBioMArray
Land.CanopyRT.create_leaf_bio
```

### Leaf optical parameters
```@docs
Land.CanopyRT.AbstractLeafOptiPara
Land.CanopyRT.LeafOptiParaArray
Land.CanopyRT.LeafOptiParaMArray
Land.CanopyRT.create_opti_par
```

### Soil
```@docs
Land.CanopyRT.SoilOpti
```

### Solar angles
```@docs
Land.CanopyRT.SolarAngles
```

### Pre-set wave length settings
```@docs
Land.CanopyRT.AbstractWLParaSet
Land.CanopyRT.WLParaSetArray
Land.CanopyRT.WLParaSetMArray
Land.CanopyRT.create_wl_para_set
```