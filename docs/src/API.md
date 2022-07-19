# ClimaCache

```@meta
CurrentModule = ClimaCache
```


## Environmental Conditions
```@docs
AirLayer
Meteorology
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
HyperspectralLeafBiophysics{FT}(wls::WaveLengthSet{FT}) where {FT<:AbstractFloat}
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
Leaf{FT}(psm::String) where {FT<:AbstractFloat}
Leaves1D
Leaves1D{FT}(psm::String) where {FT<:AbstractFloat}
Leaves2D
Leaves2D{FT}(psm::String, wls::WaveLengthSet{FT} = WaveLengthSet{FT}(); n_azi::Int = 36, n_incl::Int = 9) where {FT<:AbstractFloat}
```


## Radiation
```@docs
WaveLengthSet
WaveLengthSet{FT}(swl::Vector = WAVELENGTHS; opti::String = OPTI_2021) where {FT<:AbstractFloat}
AbstractRadiation
BroadbandRadiation
BroadbandRadiation{FT}() where {FT<:AbstractFloat}
HyperspectralRadiation
HyperspectralRadiation{FT}(wls::WaveLengthSet = WaveLengthSet{FT}(); file::String = FILE_SUN) where {FT<:AbstractFloat}
HyperspectralAbsorption
HyperspectralAbsorption{FT}(wls::WaveLengthSet = WaveLengthSet{FT}(); opti::String = OPTI_2021) where {FT<:AbstractFloat}
SunSensorGeometry
SunSensorGeometry{FT}(; haa::Number = 0, hsa::Number = 0, saa::Number = 180, sza::Number = 30, vaa::Number = 180, vza::Number = 0) where {FT<:AbstractFloat}
HyperspectralMLCanopyOpticalProperty
HyperspectralMLCanopyOpticalProperty{FT}(; n_azi::Int = 36, n_incl::Int = 9, n_layer::Int = 20, n_λ::Int = 114) where {FT<:AbstractFloat}
AbstractCanopyRadiationProfile
BroadbandSLCanopyRadiationProfile
BroadbandSLCanopyRadiationProfile{FT}(; n_incl::Int = 9) where {FT<:AbstractFloat}
HyperspectralMLCanopyRadiationProfile
HyperspectralMLCanopyRadiationProfile{FT}(; n_azi::Int = 36, n_incl::Int = 9, n_layer::Int = 20, n_λ::Int = 114) where {FT<:AbstractFloat}
AbstractLIDFAlgorithm
VerhoefLIDF
AbstractCanopy
BroadbandSLCanopy
BroadbandSLCanopy{FT}(; lai::Number = 3, θ_incl_bnds::Matrix = [collect(0:10:80) collect(10:10:90)]) where {FT<:AbstractFloat}
HyperspectralMLCanopy
HyperspectralMLCanopy{FT}(wls::WaveLengthSet{FT} = WaveLengthSet{FT}(); lai::Number = 3, n_layer::Int = 20, θ_incl_bnds::Matrix = [collect(0:10:80) collect(10:10:90)]) where {FT<:AbstractFloat}
```

## Soil
```@docs
AbstractSoilVC
BrooksCorey
VanGenuchten
VanGenuchten{FT}(name::String) where {FT<:AbstractFloat}
AbstractSoilAlbedo
BroadbandSoilAlbedo
HyperspectralSoilAlbedo
HyperspectralSoilAlbedo{FT}(wls::WaveLengthSet{FT}= WaveLengthSet{FT}()) where {FT<:AbstractFloat}
Soil
Soil{FT}(zs::Vector, area::Number, wls::WaveLengthSet{FT} = WaveLengthSet{FT}(); soil_type::String = "Loam") where {FT<:AbstractFloat}
Soil{FT}(zs::Vector, area::Number, broadband::Bool; soil_type::String = "Loam") where {FT<:AbstractFloat}
```

## SPAC
```@docs
Root
Root{FT}(; ssm::Bool = true) where {FT<:AbstractFloat}
Stem
Stem{FT}(; ssm::Bool = true) where {FT<:AbstractFloat}
AbstractSPACSystem
MonoElementSPAC
MonoElementSPAC{FT}(psm::String, zs::Vector = [-0.2,1], area::Number = 1; broadband::Bool = false, ssm::Bool = true) where {FT<:AbstractFloat}
MonoMLGrassSPAC
MonoMLGrassSPAC{FT}(psm::String, area::Number = 100, wls::WaveLengthSet{FT} = WaveLengthSet{FT}(); zs::Vector = [-0.2,0.5], zss::Vector = collect(0:-0.1:-1), zas::Vector = collect(0:0.05:1), ssm::Bool = true) where {FT<:AbstractFloat}
MonoMLPalmSPAC
MonoMLPalmSPAC{FT}(psm::String, area::Number = 100, wls::WaveLengthSet{FT} = WaveLengthSet{FT}(); zs::Vector = [-1,6,12], zss::Vector = collect(0:-0.25:-2), zas::Vector = collect(0:0.2:13), ssm::Bool = true) where {FT<:AbstractFloat}
MonoMLTreeSPAC
MonoMLTreeSPAC{FT}(psm::String, area::Number = 100, wls::WaveLengthSet{FT} = WaveLengthSet{FT}(); zs::Vector = [-1,6,12], zss::Vector = collect(0:-0.25:-2), zas::Vector = collect(0:0.5:13), ssm::Bool = true) where {FT<:AbstractFloat}
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
CP_L
CP_V
GAS_R
GRAVITY
H_PLANCK
K_BOLTZMANN
K_STEFAN
K_VON_KARMAN
LH_V0
LIGHT_SPEED
M_DRYAIR
M_H₂O
P_ATM
PRESS_TRIPLE
R_V
RT_25
T_0
T_25
T_TRIPLE
V_H₂O
YEAR_D
ρ_H₂O
ρg_MPa
```
