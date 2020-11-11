# SoilPlantAirContinuum
```@meta
CurrentModule = Land.SoilPlantAirContinuum
```

## Types
```@docs
SPACContainer1L
SPACContainer2L
SPACMono
SPACSimple
```

## Soil
```@docs
soil_moisture_swc!
soil_moisture_p!
soil_moisture_p25!
soil_moisture!
```

## Planet
```@docs
atmospheric_pressure_ratio
atmospheric_pressure
ppm_to_Pa
zenith_angle
```

## Big-leaf model
```@docs
gain_risk_map
leaf_gas_exchange_nonopt!
leaf_gas_exchange!
optimize_flows!
big_leaf_partition!
radiative_conductance
black_body_emittance
boundary_layer_conductance
leaf_temperature
leaf_temperature_sunlit
leaf_temperature_shaded
annual_profit
annual_simulation!
create_dataframe
initialize_spac_canopy!
```

## Optimal investment
```@docs
leaf_allocation!
optimize_leaf!
optimize_hs!
```
