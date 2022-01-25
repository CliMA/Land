# ClimaCache

```@meta
CurrentModule = ClimaCache
```


## Air
```@docs
AirLayer
AirLayer{FT}() where {FT<:AbstractFloat}
```


## Plant
```@docs
Leaf
Leaf{FT}(psm::String, wls::WaveLengthSet{FT} = WaveLengthSet{FT}()) where {FT<:AbstractFloat}
LeafBiophysics
LeafBiophysics{FT}(wls::WaveLengthSet{FT} = WaveLengthSet{FT}()) where {FT<:AbstractFloat}
AbstractFluorescenceModel
VanDerTolFluorescenceModel
VanDerTolFluorescenceModel{FT}(drought::Bool = false) where {FT<:AbstractFloat}
CytochromeFluorescenceModel
AbstractReactionCenter
VJPReactionCenter
VJPReactionCenter{FT}() where {FT<:AbstractFloat}
CytochromeReactionCenter
CytochromeReactionCenter{FT}() where {FT<:AbstractFloat}
```


## Radiation
```@docs
WaveLengthSet
WaveLengthSet{FT}(swl::Vector=WAVELENGTHS; opti::String=OPTI_2021) where {FT<:AbstractFloat}
```
