# Utils for the LAND model

## LandParameters module
```@docs
Land.LandParameters.EarthParameterSet
Land.LandParameters.EARTH
Land.LandParameters.AVOGADRO
Land.LandParameters.CP_D
Land.LandParameters.CP_I
Land.LandParameters.CP_L
Land.LandParameters.CP_V
Land.LandParameters.GAS_R
Land.LandParameters.GRAVITY
Land.LandParameters.H_PLANCK
Land.LandParameters.K_0
Land.LandParameters.K_25
Land.LandParameters.K_BOLTZMANN
Land.LandParameters.LH_S0
Land.LandParameters.LH_V0
Land.LandParameters.LIGHT_SPEED
Land.LandParameters.MOLMASS_DRYAIR
Land.LandParameters.MOLMASS_WATER
Land.LandParameters.PRESS_TRIPLE
Land.LandParameters.R_D
Land.LandParameters.R_V
Land.LandParameters.T_TRIPLE
Land.LandParameters.VON_KARMAN_CONST
Land.LandParameters.WATER_AIR_MRATIO
Land.LandParameters.ρ_H₂O
```

## MathTools module
```@docs
Land.MathTools.e2phot
Land.MathTools.fast∫
Land.MathTools.lower_quadratic
Land.MathTools.quadratic
Land.MathTools.volscatt
Land.MathTools.weibull_k_ratio
```

## WaterPhysics module
```@docs
Land.WaterPhysics.relative_diffusive_coefficient
Land.WaterPhysics.latent_heat_vapor
Land.WaterPhysics.saturation_vapor_pressure
Land.WaterPhysics.saturation_vapor_pressure_slope
Land.WaterPhysics.surface_tension
Land.WaterPhysics.relative_surface_tension
Land.WaterPhysics.viscosity
Land.WaterPhysics.relative_viscosity
```

## RootSolvers extension
```@docs
Land.RootSolversExtension.BisectionMethod
Land.RootSolversExtension.NewtonBisectionMethod
Land.RootSolversExtension.if_break
Land.RootSolversExtension.find_zero_ext
```
