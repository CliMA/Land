# LeafOptics

```@meta
CurrentModule = LeafOptics
```


## Leaf spectra
```@docs
leaf_spectra!
leaf_spectra!(bio::HyperspectralLeafBiophysics{FT}, wls::WaveLengthSet{FT}, lha::HyperspectralAbsorption{FT}, lwc::FT; APAR_car::Bool = true, reabsorb::Bool = true, α::FT = FT(40)) where {FT<:AbstractFloat}
leaf_spectra!(bio::HyperspectralLeafBiophysics{FT}, wls::WaveLengthSet{FT}, ρ_par::FT, ρ_nir::FT, τ_par::FT, τ_nir::FT) where {FT<:AbstractFloat}
leaf_spectra!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat}
```


## Leaf PAR, APAR, and PPAR
```@docs
leaf_PAR
```


## Leaf SIF
```@docs
leaf_SIF
```


## Utility functions
```@docs
average_transmittance
photon
photon!
energy
energy!
```
