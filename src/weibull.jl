###############################################################################
#
# Weibull functions
#
###############################################################################
"""
    weibull_k_ratio(b::FT, c::FT, p_25::FT, vis::FT)
    weibull_k_ratio(b::FT, c::FT, p::FT)

Returns the relative hydraulic conductance, given
- `b` Weibull B in `[MPa]`
- `c` Weibull C
- `p_25` Equivalent xylem pressure at 298.15 K in `[MPa]`
- `vis` Relative viscosity. If missing, vis = 1.
"""
function weibull_k_ratio(
            b::FT,
            c::FT,
            p_25::FT,
            vis::FT
            ) where {FT<:AbstractFloat}
    if p_25<0
        kr = max( FT(1e-4), exp( -1 * (-p_25/b) ^ c ) / vis );
    else
        kr = 1 / vis;
    end

    return kr
end

function weibull_k_ratio(
            b::FT,
            c::FT,
            p::FT
            ) where {FT<:AbstractFloat}
    if p<0
        kr = max( FT(1e-4), exp( -1 * (-p/b) ^ c ) );
    else
        kr = FT(1);
    end

    return kr
end
