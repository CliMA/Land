# CliMA Land Modules

- `CanopyLayers` was the canopy radiative transfer module for CliMA Land v0.1. Since CliMA Land v0.2, the core functions have been moved to new packages, and thus this module has been deprecated.
- `CanopyRadiativeTransfer` supports multiple canopy radiative transfer schemes, including multiple layered hyperspectral RT and single/multiple layered broadband RT.
- `ClimaCache` was meant to host all types and structs of CliMA Land v0.2. It has been renamed to `EmeraldNamespace` to be more accurate, and thus is deprecated till we use it for another purpose.
- `EarthSurface` contains pack of functions to compute or read constants, angles, elevations, and pressures etc at the earth surface.
- `EmeraldConstants` hosts constants used for CliMA Land model.
- `EmeraldNamespace` is designed to share information among different CliMA Land modules. All public types and structs are supposed to live here.
- `EmeraldOptics` is a package used to host the physics-based generic functions for optics.
- `LeafOptics` is a remasted package for leaf transmittance, refleclectance, and fluorescence spectra. The main functionalities are inherited from CanopyLayers.jl, which is now deprecated.
- `Photosynthesis` is designed to support freely customized photosynthesis and fluorescence models.
- `PlantHydraulics` hosts models for plant hydraulics using numerical and analytical methods.
- `SoilHydraulics` is the soil hydraulics module for Soil-Plant-Air Continuum modeling.
- `SoilPlantAirContinuum` simulates the Soil-Plant-Air Continuum with empirical and optimal stomatal models.
- `StomataModels` includes a collection of empirical and optimization stomatal models.
- `WaterPhysics` includes a collection of temperature dependencies of physical properties of water (and other traces in water).
