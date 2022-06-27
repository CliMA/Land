# LeafOptics

```@meta
CurrentModule = LeafOptics
```


## Leaf spectra
```@docs
leaf_spectra!
leaf_spectra!(bio::HyperspectralLeafBiophysics{FT}, wls::WaveLengthSet{FT}, lha::HyperspectralAbsorption{FT}; APAR_car::Bool = true, α::FT=FT(40)) where {FT<:AbstractFloat}
leaf_spectra!(bio::HyperspectralLeafBiophysics{FT}, wls::WaveLengthSet{FT}, ρ_par::FT, ρ_nir::FT, τ_par::FT, τ_nir::FT) where {FT<:AbstractFloat}
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
photon(λ::FT, E::FT) where {FT<:AbstractFloat}
photon!
photon!(λ::Vector{FT}, E::Vector{FT}, phot::Vector{FT}) where {FT<:AbstractFloat}
photon!(λ::Vector{FT}, E::Vector{FT}) where {FT<:AbstractFloat}
energy
energy(λ::FT, phot::FT) where {FT<:AbstractFloat}
energy!
energy!(λ::Vector{FT}, phot::Vector{FT}, E::Vector{FT}) where {FT<:AbstractFloat}
energy!(λ::Vector{FT}, phot::Vector{FT}) where {FT<:AbstractFloat}
```
