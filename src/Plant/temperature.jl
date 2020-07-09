"""
    radiative_conductance(T::FT) where {FT<:AbstractFloat}

Return the radiative conductance of leaf, given
- `T` leaf temperature
"""
function radiative_conductance(
    T::FT
) where {FT<:AbstractFloat}
    _A0::FT = 0.1579;
    _A1::FT = 0.0017;
    _A2::FT = 7.17E-6;
    _T ::FT = T;# - K_0;

    return _A0 + _A1*T + _A2*T^2
end

"""
    black_body_emittance(T::FT) where {FT<:AbstractFloat}

Return the energy been radiated out, given
- `T` leaf temperature
"""
function black_body_emittance(
    T::FT
) where {FT<:AbstractFloat}
    return FT(K_BOLTZMANN) * (T+FT(K_0))^4
end




"""
    boundary_layer_conductance(ws::FT, lw::FT) where {FT<:AbstractFloat}

Return the boundary layer conductance, given
- `ws` Wind speed
- `lw` Leaf width
"""
function boundary_layer_conductance(
    ws::FT,
    lw,#::FT
) where {FT<:AbstractFloat}
    _A::FT = FT(1.4 * 0.135)
    _B::FT = FT(0.72)

    return _A * sqrt( ws / (_B*lw) )
end




"""
    leaf_temperature(leaves::Leaves{FT}, envir::AirLayer{FT}, rad::FT, emol::FT, ws::FT) where {FT<:AbstractFloat}

Return leaf temperature, given
- `leaves` `Leaves` struct in `StomataModels`
- `envir` `AirLayer` struct in `Photosynthesis`
- `rad` Absorbed radiation
- `emol` Flow rate per LA
- `ws` Wind speed
"""
function leaf_temperature(
    leaves::Leaves{FT},
    envir::AirLayer{FT},
    rad::FT,
    emol::FT,
    ws::FT
) where {FT<:AbstractFloat}
    @unpack LAI,width = leaves;
    @unpack t_air = envir;

    e_emit  = FT(0.97) * black_body_emittance(t_air)
    lambda  = latent_heat_vapor(t_air) * 1000 / 18
    e_vapor = emol * lambda * LAI
    Gr      = radiative_conductance(t_air)
    GHa     = boundary_layer_conductance(ws, width)
    t_leaf  = t_air + (rad-e_vapor) / (FT(29.3)*(Gr+GHa)) / 2

    return t_leaf
end
