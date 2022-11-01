# StomtaModels API
```@meta
CurrentModule = StomataModels
```

## Empirical models
```@docs
empirical_equation
empirical_equation(sm::BallBerrySM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat}
empirical_equation(sm::BallBerrySM{FT}, leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat}
empirical_equation(sm::BallBerrySM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT<:AbstractFloat}
empirical_equation(sm::BallBerrySM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT<:AbstractFloat}
```

## Optimality models
```@docs
∂A∂E
∂R∂E
∂T∂E
∂Θ∂E
∂Θ∂E(sm::AndereggSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
∂Θ∂E(sm::AndereggSM{FT}, leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
∂Θ∂E(sm::AndereggSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
∂Θ∂E(sm::AndereggSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
∂Θₙ∂E
```

## Stomtal conductance limits
```@docs
limit_stomatal_conductance!
```

## Prognostic conductance
```@docs
∂g∂t
∂g∂t(leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
∂g∂t(leaves::Leaves1D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
∂g∂t(leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
∂g∂t(leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1), δe::FT = FT(1e-7)) where {FT<:AbstractFloat}
∂gₙ∂t
stomatal_conductance!
stomatal_conductance!(spac::MonoElementSPAC{FT}; β::FT = FT(1)) where {FT<:AbstractFloat}
stomatal_conductance!(spac::MonoElementSPAC{FT}, Δt::FT) where {FT<:AbstractFloat}
```
