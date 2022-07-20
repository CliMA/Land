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
```

### Universal Constants
```@docs
UniversalConstants
AVOGADRO
CP_D
CP_D_MOL
CP_I
CP_I_MOL
CP_L
CP_L_MOL
CP_V
CP_V_MOL
F_O₂
GAS_R
GRAVITY
H_PLANCK
K_BOLTZMANN
K_STEFAN
K_VON_KARMAN
LH_V₀
LIGHT_SPEED
M_DRYAIR
M_H₂O
P_ATM
PRESS_TRIPLE
R_V
RT₂₅
T₀
T₂₅
T_TRIPLE
V_H₂O
YEAR_D
Λ_THERMAL_H₂O
ρ_H₂O
ρg_MPa
```
