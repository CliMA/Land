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
Leaf{FT}(psm::String, wls::WaveLengthSet{FT} = WaveLengthSet{FT}(); colimit::Bool = false) where {FT<:AbstractFloat}
LeafBiophysics
LeafBiophysics{FT}(wls::WaveLengthSet{FT} = WaveLengthSet{FT}()) where {FT<:AbstractFloat}
VanDerTolFluorescenceModel
VanDerTolFluorescenceModel{FT}(drought::Bool = false) where {FT<:AbstractFloat}
AbstractReactionCenter
VJPReactionCenter
VJPReactionCenter{FT}() where {FT<:AbstractFloat}
CytochromeReactionCenter
CytochromeReactionCenter{FT}() where {FT<:AbstractFloat}
AbstractPhotosynthesisModel
C3CytochromeModel
C3CytochromeModel{FT}(; v_cmax25::Number = 50, r_d25::Number = 0.75, colimit::Bool = false) where {FT<:AbstractFloat}
C3VJPModel
C3VJPModel{FT}(; v_cmax25::Number = 50, j_max25::Number = 83.5, r_d25::Number = 0.75, colimit::Bool = false) where {FT<:AbstractFloat}
C4VJPModel
C4VJPModel{FT}(; v_cmax25::Number = 50, v_pmax25::Number = 50, r_d25::Number = 0.75, colimit::Bool = false) where {FT<:AbstractFloat}
AbstractPhotosynthesisMode
GCO₂Mode
PCO₂Mode
AbstractColimit
MinimumColimit
QuadraticColimit
SerialColimit
AbstractTemperatureDependency
Arrhenius
ArrheniusPeak
Q10
```


## Radiation
```@docs
WaveLengthSet
WaveLengthSet{FT}(swl::Vector=WAVELENGTHS; opti::String=OPTI_2021) where {FT<:AbstractFloat}
```
