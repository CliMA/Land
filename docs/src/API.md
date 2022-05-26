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

### Hydraulics
```@docs
AbstractPVCurve
LinearPVCurve
LinearPVCurve{FT}() where {FT<:AbstractFloat}
SegmentedPVCurve
SegmentedPVCurve{FT}() where {FT<:AbstractFloat}
AbstractHydraulicSystem
AbstractXylemVC
LogisticVC
PowerVC
WeibullVC
ComplexVC
AbstractSteadyStateMode
SteadyStateMode
NonSteadyStateMode
LeafHydraulics
LeafHydraulics{FT}(N::Int = 5; area::Number = 1500, k_ox::Number = 100, k_sla::Number = 0.04, v_max::Number = 20) where {FT<:AbstractFloat}
RootHydraulics
RootHydraulics{FT}(N::Int = 5; area::Number = 1, k_x::Number = 25, Δh = 1) where {FT<:AbstractFloat}
```

### Leaf Level
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
WaveLengthSet{FT}(swl::Vector = WAVELENGTHS; opti::String = OPTI_2021) where {FT<:AbstractFloat}
HyperspectralRadiation
HyperspectralRadiation{FT}(wls::WaveLengthSet = WaveLengthSet{FT}(); file::String = FILE_SUN) where {FT<:AbstractFloat}
HyperspectralAbsorption
HyperspectralAbsorption{FT}(wls::WaveLengthSet = WaveLengthSet{FT}(); opti::String = OPTI_2021) where {FT<:AbstractFloat}
```

## Soil
```@docs
AbstractSoilVC
BrooksCorey
VanGenuchten
VanGenuchten{FT}(name::String, α::Number, n::Number, θ_sat::Number, θ_res::Number) where {FT<:AbstractFloat}
VanGenuchten{FT}(name::String) where {FT<:AbstractFloat}
```
