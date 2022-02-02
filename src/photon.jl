# constants
const FAC = 1e-9 / (H_PLANCK() * LIGHT_SPEED() * AVOGADRO());


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Oct-22: rename the function to photon
#     2022-Feb-02: fix documentation
#
#######################################################################################################################################################################################################
"""
This function convert energy to photons based the wavelength. Supported methods are

$(METHODLIST)

"""
function photon end


#######################################################################################################################################################################################################
#
# Changes made to this method
# General
#     2021-Oct-22: add a method to convert direct from number to number
#     2022-Feb-02: fix documentation
#
#######################################################################################################################################################################################################
"""

    photon(λ::FT, E::FT) where {FT<:AbstractFloat}

Return the number of moles of photons, given
- `λ` Wave length in `[nm]`, converted to `[m]` by FAC
- `E` Joules of energy

---
# Examples
```julia
phot = photon(500.0, 0.1);
phots = photon.([500.0, 400.0], [0.1, 0.1]);
```
"""
photon(λ::FT, E::FT) where {FT<:AbstractFloat} = E * λ * FT(FAC);


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Oct-22: define function to convert photon back to energy
#     2022-Feb-02: fix documentation
#
#######################################################################################################################################################################################################
"""
This function convert photons to energy based the wavelength. Supported methods are

$(METHODLIST)

"""
function energy end


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Oct-22: add a method to convert photon back to energy
#     2022-Feb-02: fix documentation
#
#######################################################################################################################################################################################################
"""

    energy(λ::FT, phot::FT) where {FT<:AbstractFloat}

Return the energy, given
- `λ` Wave length in `[nm]`, converted to `[m]` by FAC
- `phot` Number of moles of photon

---
# Examples
```julia
e = energy(500.0, 0.1);
es = energy.([500.0, 400.0], [0.1, 0.1]);
```
"""
energy(λ::FT, phot::FT) where {FT<:AbstractFloat} = phot / (λ * FT(FAC));
