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
