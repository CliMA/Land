###############################################################################
#
# psofunction
#
###############################################################################
"""
    psofunction(K::FT,
                k::FT,
                Ω::FT,
                LAI::FT,
                q::FT,
                dso::FT,
                xl::FT
    ) where {FT<:AbstractFloat}

# TODO explain the variables
Return the probability of observing a sunlit leaf at depth `xl` (`pso`, see eq
    31 in vdT 2009), given
- `xl` Leaf depth in the canopy
"""
function psofunction(
            K::FT,
            k::FT,
            Ω::FT,
            LAI::FT,
            q::FT,
            dso::FT,
            xl::FT
) where {FT<:AbstractFloat}
    # [nl+1]  factor for correlation of Ps and Po
    if dso != 0
        alf = dso / q * 2 / (k + K);

        return exp( (K+k)*Ω*LAI*xl +
                    sqrt(K*k)*Ω*LAI / alf * (1 - exp(xl*alf)) )
    else
        return exp( (K + k -sqrt(K*k)) * Ω * LAI * xl )
    end
end
