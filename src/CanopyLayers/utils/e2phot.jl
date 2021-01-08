###############################################################################
#
# Convert energy to photons
#
###############################################################################
_FAC(FT) = FT(1e-9) / (H_PLANCK(FT) * LIGHT_SPEED(FT) * AVOGADRO(FT))

"""
    e2phot(λ::Array{FT}, E::Array{FT}) where {FT<:AbstractFloat}

Calculates the number of moles of photons, given
- `λ` An array of wave length in `[nm]`, converted to `[m]` by _FAC
- `E` Joules of energy
"""
function e2phot(
            λ::Array{FT},
            E::Array{FT}
) where {FT<:AbstractFloat}
    return (E .* λ) .* _FAC(FT)
end




"""
    e2phot!(
                λ::Array{FT,1},
                E::Array{FT,1},
                cache::Array{FT,1}
    ) where {FT<:AbstractFloat}

Calculates the number of moles of photons, given
- `λ` An array of wave length in `[nm]`, converted to `[m]` by _FAC
- `E` Joules of energy
- `cache` Cache to avoid memory allocations
"""
function e2phot!(
            λ::Array{FT,1},
            E::Array{FT,1},
            cache::Array{FT,1}
) where {FT<:AbstractFloat}
    cache .= (E .* λ) .* _FAC(FT);

    return nothing
end
