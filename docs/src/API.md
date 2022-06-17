# CanopyRadiativeTransfer

```@meta
CurrentModule = CanopyRadiativeTransfer
```

## Leaf Inclination Angle
```@docs
lidf_cdf
lidf_cdf(lidf::VerhoefLIDF{FT}, θ::FT) where {FT<:AbstractFloat}
inclination_angles!
inclination_angles!(can::Union{BroadbandSLCanopy{FT}, HyperspectralMLCanopy{FT}}, lidf::VerhoefLIDF{FT}) where {FT<:AbstractFloat}
inclination_angles!(can::Union{BroadbandSLCanopy{FT}, HyperspectralMLCanopy{FT}}, lidf::VerhoefLIDF{FT}, a::FT, b::FT) where {FT<:AbstractFloat}
```

## Clumping Index
```@docs
clumping_index!
```

## Hyperspectral Canopy RT
```@docs
soil_albedo!
soil_albedo!(can::HyperspectralMLCanopy{FT}, soil::Soil{FT}, albedo::BroadbandSoilAlbedo{FT}; clm::Bool = false) where {FT<:AbstractFloat}
soil_albedo!(can::HyperspectralMLCanopy{FT}, soil::Soil{FT}, albedo::HyperspectralSoilAlbedo{FT}; clm::Bool = false) where {FT<:AbstractFloat}
soil_albedo!(can::HyperspectralMLCanopy{FT}, soil::Soil{FT}) where {FT<:AbstractFloat}
extinction_coefficient
extinction_coefficient(sza::FT, lia::FT) where {FT<:AbstractFloat}
extinction_coefficient(lia::FT) where {FT<:AbstractFloat}
extinction_coefficient(sza::FT, vza::FT, raa::FT, lia::FT) where {FT<:AbstractFloat}
extinction_scattering_coefficients!
extinction_scattering_coefficients!(can::BroadbandSLCanopy{FT}, angles::SunSensorGeometry{FT}) where {FT<:AbstractFloat}
extinction_scattering_coefficients!(can::HyperspectralMLCanopy{FT}, angles::SunSensorGeometry{FT}) where {FT<:AbstractFloat}
canopy_optical_properties!
canopy_optical_properties!(can::HyperspectralMLCanopy{FT}, angles::SunSensorGeometry{FT}) where {FT<:AbstractFloat}
canopy_optical_properties!(can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaf{FT}}, soil::Soil{FT}) where {FT<:AbstractFloat}
canopy_radiation!
canopy_radiation!(can::BroadbandSLCanopy{FT}, leaf::Leaf{FT}, rad::BroadbandRadiation{FT}) where {FT<:AbstractFloat}
canopy_radiation!(can::HyperspectralMLCanopy{FT}, albedo::BroadbandSoilAlbedo{FT}) where {FT<:AbstractFloat}
canopy_radiation!(can::HyperspectralMLCanopy{FT}, albedo::HyperspectralSoilAlbedo{FT}) where {FT<:AbstractFloat}
canopy_radiation!(can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaf{FT}}, rad::HyperspectralRadiation{FT}, soil::Soil{FT}) where {FT<:AbstractFloat}
canopy_radiation!(can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaf{FT}}, rad::FT, soil::Soil{FT}) where {FT<:AbstractFloat}
canopy_fluorescence!
canopy_fluorescence!(can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaf{FT}}; ϕ_photon::Bool = true) where {FT<:AbstractFloat}
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
