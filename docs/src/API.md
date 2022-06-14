# CanopyRadiativeTransfer

```@meta
CurrentModule = CanopyRadiativeTransfer
```

## Leaf Inclination Angle
```@docs
lidf_cdf
lidf_cdf(lidf::VerhoefLIDF{FT}, θ::FT) where {FT<:AbstractFloat}
inclination_angles!
inclination_angles!(can::HyperspectralMLCanopy{FT}, lidf::VerhoefLIDF{FT}) where {FT<:AbstractFloat}
inclination_angles!(can::HyperspectralMLCanopy{FT}, lidf::VerhoefLIDF{FT}, a::FT, b::FT) where {FT<:AbstractFloat}
```

## Clumping Index
```@docs
clumping_index!
```

## Hyperspectral Canopy RT
```@docs
extinction_scattering_coefficients
extinction_scattering_coefficients!
canopy_optical_properties!
canopy_optical_properties!(can::HyperspectralMLCanopy{FT}, angles::SunSensorGeometry{FT}) where {FT<:AbstractFloat}
canopy_optical_properties!(can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaf{FT}}, soil::Soil{FT}) where {FT<:AbstractFloat}
canopy_radiation!
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
```
