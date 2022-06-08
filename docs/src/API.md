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
AbstractFlowProfile
NonSteadyStateFlow
NonSteadyStateFlow{FT}(N::Int, isleaf::Bool = true) where {FT<:AbstractFloat}
SteadyStateFlow
LeafHydraulics
LeafHydraulics{FT}(N::Int = 5; area::Number = 1500, k_ox::Number = 100, k_sla::Number = 0.04, v_max::Number = 20, ssm::Bool = true) where {FT<:AbstractFloat}
RootHydraulics
RootHydraulics{FT}(N::Int = 5; area::Number = 1, k_x::Number = 25, Δh::Number = 1, Δl::Number = 1, ssm::Bool = true) where {FT<:AbstractFloat}
StemHydraulics
StemHydraulics{FT}(N::Int = 5; area::Number = 1, k_x::Number = 25, Δh::Number = 1, Δl::Number = 1, ssm::Bool = true) where {FT<:AbstractFloat}
```

### Leaf Level
```@docs
Leaf
Leaf{FT}(psm::String, wls::WaveLengthSet{FT} = WaveLengthSet{FT}(); colimit::Bool = false, ssm::Bool = true) where {FT<:AbstractFloat}
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
SunSensorGeometry
SunSensorGeometry{FT}(; haa::Number = 0, hsa::Number = 0, saa::Number = 180, sza::Number = 30, vaa::Number = 180, vza::Number = 0) where {FT<:AbstractFloat}
CanopyOpticalProperty
CanopyOpticalProperty{FT}(; n_azi::Int = 36, n_incl::Int = 9, n_layer::Int = 20, n_λ::Int = 114) where {FT<:AbstractFloat}
AbstractLIDFAlgorithm
VerhoefLIDF
AbstractCanopyStructure
HyperspectralMLCanopy
HyperspectralMLCanopy{FT}(; lai::Number = 3, n_layer::Int = 20, n_λ::Int = 114, θ_incl_bnds::Matrix = [collect(0:10:80) collect(10:10:90)]) where {FT<:AbstractFloat}
```

## Soil
```@docs
AbstractSoilVC
BrooksCorey
VanGenuchten
VanGenuchten{FT}(name::String, α::Number, n::Number, θ_sat::Number, θ_res::Number) where {FT<:AbstractFloat}
VanGenuchten{FT}(name::String) where {FT<:AbstractFloat}
Soil
Soil{FT}(zs::Vector{FT}; soil_type::String = "Loam", n_λ::Int = 114) where {FT<:AbstractFloat}
```

## SPAC
```@docs
Root
Root{FT}(; ssm::Bool = true) where {FT<:AbstractFloat}
Stem
Stem{FT}(; ssm::Bool = true) where {FT<:AbstractFloat}
AbstractSPACSystem
MonoElementSPAC
MonoElementSPAC{FT}(psm::String) where {FT<:AbstractFloat}
MonoGrassSPAC
MonoGrassSPAC{FT}(psm::String; zr::Number = -0.2, zc::Number = 0.5, zss::Vector = collect(0:-0.1:-1), zas::Vector = collect(0:0.05:1), ssm::Bool = true) where {FT<:AbstractFloat}
MonoPalmSPAC
MonoPalmSPAC{FT}(psm::String; zr::Number = -1, zt::Number = 10, zc::Number = 12, zss::Vector = collect(0:-0.25:-2), zas::Vector = collect(0:0.2:13), ssm::Bool = true) where {FT<:AbstractFloat}
MonoTreeSPAC
MonoTreeSPAC{FT}(psm::String; zr::Number = -1, zt::Number = 10, zc::Number = 12, zss::Vector = collect(0:-0.25:-2), zas::Vector = collect(0:0.2:13), ssm::Bool = true) where {FT<:AbstractFloat}
```
