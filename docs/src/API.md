# Photosynthesis

```@meta
CurrentModule = Photosynthesis
```


## Photosynthesis model
```@docs
leaf_photosynthesis!
leaf_photosynthesis!(leaf::Leaf{FT}, air::AirLayer{FT}, mode::PCO₂Mode, p_i::FT = leaf.p_CO₂_i) where {FT<:AbstractFloat}
leaf_photosynthesis!(leaf::Leaf{FT}, air::AirLayer{FT}, mode::GCO₂Mode, g_lc::FT = leaf.g_CO₂) where {FT<:AbstractFloat}
```


## Temperature dependency
```@docs
photosystem_temperature_dependence!
photosystem_temperature_dependence!(psm::C3VJPModel{FT}, air::AirLayer{FT}, t::FT) where {FT<:AbstractFloat}
photosystem_temperature_dependence!(psm::C3CytochromeModel{FT}, air::AirLayer{FT}, t::FT) where {FT<:AbstractFloat}
photosystem_temperature_dependence!(psm::C4VJPModel{FT}, air::AirLayer{FT}, t::FT) where {FT<:AbstractFloat}
temperature_correction
temperature_correction(td::Arrhenius{FT}, t::FT; t_ref::FT = td.T_REF) where {FT<:AbstractFloat}
temperature_correction(td::ArrheniusPeak{FT}, t::FT; t_ref::FT = td.T_REF) where {FT<:AbstractFloat}
temperature_correction(td::Q10{FT}, t::FT; t_ref::FT = td.T_REF) where {FT<:AbstractFloat}
temperature_corrected_value
```


## Electron transport
```@docs
photosystem_electron_transport!
photosystem_electron_transport!(psm::C3VJPModel{FT}, rc::VJPReactionCenter{FT}, apar::FT, p_i::FT) where {FT<:AbstractFloat}
photosystem_electron_transport!(psm::C3CytochromeModel{FT}, rc::CytochromeReactionCenter{FT}, apar::FT, p_i::FT) where {FT<:AbstractFloat}
photosystem_electron_transport!(psm::C4VJPModel{FT}, rc::VJPReactionCenter{FT}, apar::FT, p_i::FT) where {FT<:AbstractFloat}
```


## Photosynthetic rates
```@docs
rubisco_limited_rate!
rubisco_limited_rate!(psm::Union{C3CytochromeModel{FT},C3VJPModel{FT}}, p_i::FT) where {FT<:AbstractFloat}
rubisco_limited_rate!(psm::C4VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat}
rubisco_limited_rate!(psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat}
rubisco_limited_rate!(psm::C4VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat}
light_limited_rate!
light_limited_rate!(psm::C3VJPModel{FT}) where {FT<:AbstractFloat}
light_limited_rate!(psm::Union{C3CytochromeModel{FT}, C4VJPModel{FT}}) where {FT<:AbstractFloat}
light_limited_rate!(psm::C3CytochromeModel{FT}, rc::CytochromeReactionCenter{FT}, air::AirLayer{FT}, apar::FT, g_lc::FT) where {FT<:AbstractFloat}
light_limited_rate!(psm::C3VJPModel{FT}, rc::VJPReactionCenter{FT}, air::AirLayer{FT}, apar::FT, g_lc::FT) where {FT<:AbstractFloat}
light_limited_rate!(psm::C4VJPModel{FT}, rc::VJPReactionCenter{FT}, air::AirLayer{FT}, apar::FT, g_lc::FT) where {FT<:AbstractFloat}
product_limited_rate!
product_limited_rate!(psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}}, p_i::FT) where {FT<:AbstractFloat}
product_limited_rate!(psm::C4VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat}
product_limited_rate!(psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat}
product_limited_rate!(psm::C4VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat}
```


## Colimitation
```@docs
colimit_photosynthesis!
colimit_photosynthesis!(psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}, C4VJPModel{FT}}) where {FT<:AbstractFloat}
colimit_photosynthesis!(psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}, C4VJPModel{FT}}, colim_cj::Union{MinimumColimit{FT}, QuadraticColimit{FT}},colim_ip::Union{MinimumColimit{FT}, QuadraticColimit{FT}}) where {FT<:AbstractFloat}
colimited_rate
colimited_rate(a_1::FT, a_2::FT, colim::MinimumColimit{FT}) where {FT<:AbstractFloat}
colimited_rate(a_1::FT, a_2::FT, colim::QuadraticColimit{FT}) where {FT<:AbstractFloat}
```


## Coefficients and fluorescence
```@docs
photosystem_coefficients!
photosystem_coefficients!(psm::Union{C3VJPModel{FT}, C4VJPModel{FT}}, rc::VJPReactionCenter{FT}, apar::FT) where {FT<:AbstractFloat}
photosystem_coefficients!(psm::C3CytochromeModel{FT}, rc::CytochromeReactionCenter{FT}, apar::FT) where {FT<:AbstractFloat}
```
