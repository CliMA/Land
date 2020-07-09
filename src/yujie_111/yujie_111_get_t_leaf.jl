# assume air emissivity is 0.3 and that for leaf is 0.97
function Yujie111GetLeafTem(node::Yujie111{FT}, tair, rad, emol, shading, wind) where {FT<:AbstractFloat}
    e_emit  = (0.97-0.3) * black_body_emittance(tair)
    lambda  = -42.9143*tair + 45064.3
    e_vapor = emol * lambda / node.gaba
    Gr      = radiative_conductance(tair)
    GHa     = boundary_layer_conductance(wind, node.width)
    if shading==false
        t_leaf  = tair + 0.5(rad-e_emit-e_vapor) / (29.3*(Gr+GHa))
    else
        t_leaf  = tair + 0.5(rad-e_vapor) / (29.3*(Gr+GHa))
    end
    return t_leaf
end
