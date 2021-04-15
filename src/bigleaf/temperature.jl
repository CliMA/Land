###############################################################################
#
# Radiative conductance
#
###############################################################################
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
    _T ::FT = T - T_0(FT);

    return _A0 + _A1*T + _A2*T^2
end








###############################################################################
#
# Black body energy emittance
#
###############################################################################
"""
    black_body_emittance(T::FT) where {FT<:AbstractFloat}

Return the energy been radiated out, given
- `T` leaf temperature
"""
function black_body_emittance(
            T::FT
) where {FT<:AbstractFloat}
    return K_STEFAN(FT) * T^4
end








###############################################################################
#
# Leaf boundary layer conductance
#
###############################################################################
"""
    boundary_layer_conductance(wind::FT, width::FT) where {FT<:AbstractFloat}

Return the boundary layer conductance, given
- `wind` Wind speed
- `width` Leaf width
"""
function boundary_layer_conductance(
            wind::FT,
            width::FT
) where {FT<:AbstractFloat}
    _A::FT = FT(1.4 * 0.135)
    _B::FT = FT(0.72)

    return _A * sqrt( wind / (_B*width) )
end








###############################################################################
#
# Leaf temperature
#
###############################################################################
"""
    leaf_temperature(
                node::SPACSimple{FT},
                rad::FT,
                e_rad::FT,
                epla::FT
    ) where {FT<:AbstractFloat}
    leaf_temperature(
                node::SPACSimple{FT},
                rad::FT,
                flow::FT
    ) where {FT<:AbstractFloat}

Return leaf temperature, given
- `node` [`SPACSimple`](@ref) type struct
- `rad` Absorbed solar radiation per leaf area
- `e_rad` Emitted absorbed radiation per leaf area
- `epla` Flow rate per leaf area
- `flow` Total flow rate per basal area
"""
function leaf_temperature(
            node::SPACSimple{FT},
            rad::FT,
            e_rad::FT,
            epla::FT
) where {FT<:AbstractFloat}
    @unpack t_air,wind = node.envir;

    lambda = latent_heat_vapor(t_air) * M_H₂O(FT);
    Gr     = radiative_conductance(t_air);
    GHa    = boundary_layer_conductance(wind, node.width);
    e_lat  = lambda * epla;

    t_leaf = t_air + (rad-e_rad-e_lat) / (FT(29.3)*(Gr+GHa)) / 2;

    return t_leaf
end




function leaf_temperature(
            node::SPACSimple{FT},
            rad::FT,
            flow::FT
) where {FT<:AbstractFloat}
    lai = node.laba / node.gaba;

    e_rad = FT(0.97) * black_body_emittance(node.envir.t_air) / max(1, lai);
    epla  = flow / node.laba;

    return leaf_temperature(node, rad, e_rad, epla)
end




"""
    leaf_temperature_sunlit(
                node::SPACSimple{FT},
                rad::FT,
                f_sl::FT
    ) where {FT<:AbstractFloat}

Return leaf temperature, given
- `node` [`SPACSimple`](@ref) type struct
- `rad` Absorbed solar radiation per leaf area
- `f_sl` Total flow rate per basal area into sunlit leaves
"""
function leaf_temperature_sunlit(
            node::SPACSimple{FT},
            rad::FT,
            f_sl::FT
) where {FT<:AbstractFloat}
    e_rad  = FT(0.97) * black_body_emittance(node.envir.t_air) /
             max(1, (node.container2L).lai_sl);
    epla   = f_sl / (node.container2L).la_sl;

    return leaf_temperature(node, rad, e_rad, epla)
end




"""
    leaf_temperature_shaded(
                node::SPACSimple{FT},
                rad::FT,
                f_sh::FT
    ) where {FT<:AbstractFloat}

Return leaf temperature, given
- `node` [`SPACSimple`](@ref) type struct
- `rad` Absorbed solar radiation per leaf area
- `f_sh` Total flow rate per basal area into shaded leaves
"""
function leaf_temperature_shaded(
            node::SPACSimple{FT},
            rad::FT,
            f_sh::FT
) where {FT<:AbstractFloat}
    epla = f_sh / (node.container2L).la_sh;

    return leaf_temperature(node, rad, FT(0), epla)
end
