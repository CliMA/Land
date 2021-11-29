# ClimaCache
```@meta
CurrentModule = ClimaCache
```


## Plant
```@docs
HyperspectralAbsorption
HyperspectralAbsorption{FT}(wls::WaveLengthSet) where {FT<:AbstractFloat}
LeafBiophysics
LeafBiophysics{FT}(wls::WaveLengthSet) where {FT<:AbstractFloat}
AbstractPhotosynthesisSystem
C₃VJPSystem
C₃VJPSystem{FT}() where {FT<:AbstractFloat}
```


## Radiation
```@docs
HyperspectralRadiation
HyperspectralRadiation{FT}(wls::WaveLengthSet) where {FT<:AbstractFloat}
WaveLengthSet
WaveLengthSet{FT}(swl::Vector) where {FT<:AbstractFloat}
```


## Soil
```@docs
AbstractSoilVC
BrooksCorey
VanGenuchten
VanGenuchten{FT}(name::String, α::Number, n::Number, θ_sat::Number, θ_res::Number) where {FT<:AbstractFloat}
VanGenuchten{FT}(name::String) where {FT<:AbstractFloat}
```
