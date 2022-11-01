# ClimaCache

```@meta
CurrentModule = ClimaCache
```


## Environmental Conditions
```@docs
AirLayer
Meteorology
AbstractSoilVC
BrooksCorey
VanGenuchten
VanGenuchten{FT}(name::String) where {FT<:AbstractFloat}
AbstractSoilAlbedo
BroadbandSoilAlbedo
HyperspectralSoilAlbedo
SoilLayer
Soil
```


## Plant

### Hydraulics
```@docs
AbstractPVCurve
LinearPVCurve
SegmentedPVCurve
AbstractHydraulicSystem
AbstractXylemVC
LogisticVC
PowerVC
WeibullVC
ComplexVC
AbstractFlowProfile
NonSteadyStateFlow
SteadyStateFlow
LeafHydraulics
RootHydraulics
StemHydraulics
```

### Leaf Level
```@docs
AbstractLeafBiophysics
BroadbandLeafBiophysics
HyperspectralLeafBiophysics
VanDerTolFluorescenceModel
VDTModelAll
VDTModelDrought
AbstractReactionCenter
CytochromeReactionCenter
VJPReactionCenter
AbstractPhotosynthesisModel
C3CytochromeModel
C3VJPModel
C4VJPModel
AbstractPhotosynthesisMode
GCO₂Mode
PCO₂Mode
AbstractTemperatureDependency
Arrhenius
ArrheniusPeak
Q10
Q10Peak
AbstractLeaf
Leaf
Leaves1D
Leaves2D
```

### Stomatal Models
```@docs
AbstractBetaParameter
BetaParameterG1
BetaParameterKleaf
BetaParameterKsoil
BetaParameterPleaf
BetaParameterPsoil
BetaParameterVcmax
BetaParameterΘ
BetaFunction
AbstractStomataModel
AndereggSM
BallBerrySM
EllerSM
GentineSM
LeuningSM
MedlynSM
SperrySM
WangSM
Wang2SM
```

## Radiation
```@docs
WaveLengthSet
AbstractRadiation
BroadbandRadiation
HyperspectralRadiation
HyperspectralAbsorption
SunSensorGeometry
HyperspectralMLCanopyOpticalProperty
AbstractCanopyRadiationProfile
BroadbandSLCanopyRadiationProfile
HyperspectralMLCanopyRadiationProfile
AbstractLIDFAlgorithm
VerhoefLIDF
AbstractCanopy
BroadbandSLCanopy
HyperspectralMLCanopy
```

## SPAC
```@docs
Root
Stem
AbstractSPACSystem
MonoElementSPAC
MonoMLGrassSPAC
MonoMLPalmSPAC
MonoMLTreeSPAC
```


## Utils

### Colimitation Methods
```@docs
AbstractColimit
MinimumColimit
QuadraticColimit
SerialColimit
SquareColimit
```
