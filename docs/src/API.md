# LeafOptics

```@meta
CurrentModule = LeafOptics
```


## Leaf spectra
```@docs
leaf_spectra!
leaf_spectra!(bio::LeafBiophysics{FT}, wls::WaveLengthSet{FT}, lha::HyperspectralAbsorption{FT}; APAR_car::Bool = true, α::FT=FT(40)) where {FT<:AbstractFloat}
leaf_spectra!(bio::LeafBiophysics{FT}, wls::WaveLengthSet{FT}, ρ_par::FT, ρ_nir::FT, τ_par::FT, τ_nir::FT) where {FT<:AbstractFloat}
```


## Utility functions
```@docs
average_transmittance(α::FT, nr::FT) where {FT<:AbstractFloat}
photon
photon(λ::FT, E::FT) where {FT<:AbstractFloat}
energy
energy(λ::FT, phot::FT) where {FT<:AbstractFloat}
```
