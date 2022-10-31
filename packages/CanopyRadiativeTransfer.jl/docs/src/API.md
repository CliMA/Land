# CanopyRadiativeTransfer

```@meta
CurrentModule = CanopyRadiativeTransfer
```

## Leaf Inclination Angle
```@docs
lidf_cdf
inclination_angles!
```

## Clumping Index
```@docs
clumping_index!
```

## Hyperspectral Canopy RT
```@docs
soil_albedo!
extinction_coefficient
extinction_coefficient(sza::FT, lia::FT) where {FT<:AbstractFloat}
extinction_coefficient(lia::FT) where {FT<:AbstractFloat}
extinction_coefficient(sza::FT, vza::FT, raa::FT, lia::FT) where {FT<:AbstractFloat}
extinction_scattering_coefficients!
extinction_scattering_coefficients!(can::BroadbandSLCanopy{FT}, angles::SunSensorGeometry{FT}) where {FT<:AbstractFloat}
extinction_scattering_coefficients!(can::HyperspectralMLCanopy{FT}, angles::SunSensorGeometry{FT}) where {FT<:AbstractFloat}
canopy_optical_properties!
canopy_optical_properties!(can::HyperspectralMLCanopy{FT}, angles::SunSensorGeometry{FT}) where {FT<:AbstractFloat}
canopy_optical_properties!(can::HyperspectralMLCanopy{FT}, albedo::BroadbandSoilAlbedo{FT}) where {FT<:AbstractFloat}
canopy_optical_properties!(can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaves2D{FT}}, soil::Soil{FT}) where {FT<:AbstractFloat}
canopy_radiation!
canopy_radiation!(can::BroadbandSLCanopy{FT}, leaf::Leaves1D{FT}, rad::BroadbandRadiation{FT}, soil::Soil{FT}) where {FT<:AbstractFloat}
canopy_radiation!(can::HyperspectralMLCanopy{FT}, albedo::BroadbandSoilAlbedo{FT}) where {FT<:AbstractFloat}
canopy_radiation!(can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaves2D{FT}}, rad::HyperspectralRadiation{FT}, soil::Soil{FT}) where {FT<:AbstractFloat}
canopy_radiation!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat}
canopy_fluorescence!
```

## Remote Sensing Applications
```@docs
read_spectrum
read_spectrum(x::Vector{FT}, y::Vector{FT}, target::FT) where {FT<:AbstractFloat}
read_spectrum(x::Vector{FT}, y::Vector{FT}, x₁::FT, x₂::FT; steps::Int = 2) where {FT<:AbstractFloat}
MODIS_EVI
MODIS_EVI2
MODIS_LSWI
MODIS_NDVI
MODIS_NIRv
OCO2_SIF759
OCO2_SIF770
OCO3_SIF759
OCO3_SIF770
TROPOMI_SIF683
TROPOMI_SIF740
TROPOMI_SIF747
TROPOMI_SIF751
```
