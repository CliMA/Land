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
