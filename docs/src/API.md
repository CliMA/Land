# Photosynthesis

```@meta
CurrentModule = Photosynthesis
```


## Photosynthesis model
```@docs
leaf_photosynthesis!
leaf_photosynthesis!(lf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}, g_lc::FT, ppar::FT) where {FT<:AbstractFloat}
leaf_photosynthesis!(lf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT<:AbstractFloat}
leaf_photosynthesis!(spac::MonoElementSPAC{FT}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT<:AbstractFloat}
```


## Temperature dependency
```@docs
temperature_correction
temperature_corrected_value
photosystem_temperature_dependence!
∂R∂T
```


## Electron transport
```@docs
photosystem_electron_transport!
```


## Photosynthetic rates
```@docs
rubisco_limited_rate!
rubisco_limited_rate!(psm::Union{C3CytochromeModel{FT},C3VJPModel{FT}}, p_i::FT; β::FT = FT(1)) where {FT<:AbstractFloat}
rubisco_limited_rate!(psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT<:AbstractFloat}
light_limited_rate!
light_limited_rate!(psm::Union{C3CytochromeModel{FT}, C4VJPModel{FT}}) where {FT<:AbstractFloat}
light_limited_rate!(psm::C3CytochromeModel{FT}, rc::CytochromeReactionCenter{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT<:AbstractFloat}
product_limited_rate!
product_limited_rate!(psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}}, p_i::FT; β::FT = FT(1)) where {FT<:AbstractFloat}
product_limited_rate!(psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT<:AbstractFloat}
```


## Colimitation
```@docs
colimit_photosynthesis!
colimited_rate
```


## Coefficients and fluorescence
```@docs
photosystem_coefficients!
```
