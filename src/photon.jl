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
#     2022-Jun-13: add function
#     2022-Jun-15: fix documentation
#
#######################################################################################################################################################################################################
"""
This function convert energy to photons based the wavelength and save it to provided variable (1D+). Supported methods are

$(METHODLIST)

"""
function photon! end


#######################################################################################################################################################################################################
#
# Changes made to this method
# General
#     2021-Jun-13: add method to save to provided 3rd variable
#
#######################################################################################################################################################################################################
"""

    photon!(λ::Vector{FT}, E::Vector{FT}, phot::Vector{FT}) where {FT<:AbstractFloat}

Compute and save the number of moles of photons, given
- `λ` Wave length in `[nm]`, converted to `[m]` by FAC
- `E` Joules of energy
- `phot` Mole of photons (variable to save)

---
# Examples
```julia
phots = zeros(2);
photon!.([500.0, 400.0], [0.1, 0.1], phots);
```
"""
photon!(λ::Vector{FT}, E::Vector{FT}, phot::Vector{FT}) where {FT<:AbstractFloat} = (
    phot .= photon.(λ, E);

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes made to this method
# General
#     2021-Jun-13: add method to save to provided 2rd variable
#
#######################################################################################################################################################################################################
"""

    photon!(λ::Vector{FT}, E::Vector{FT}) where {FT<:AbstractFloat}

Compute and save the number of moles of photons, given
- `λ` Wave length in `[nm]`, converted to `[m]` by FAC
- `E` Joules of energy, will be converted to moles of photons

---
# Examples
```julia
Es = rand(2);
photon!.([500.0, 400.0], Es);
```
"""
photon!(λ::Vector{FT}, E::Vector{FT}) where {FT<:AbstractFloat} = (
    E .*= λ .* FT(FAC);

    return nothing
);


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


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2022-Jun-13: add function
#     2022-Jun-15: fix documentation
#
#######################################################################################################################################################################################################
"""
This function convert photons to energy based the wavelength and save it to provided variable (1D+). Supported methods are

$(METHODLIST)

"""
function energy! end


#######################################################################################################################################################################################################
#
# Changes made to this method
# General
#     2021-Jun-13: add method to save to provided 3rd variable
#
#######################################################################################################################################################################################################
"""

    energy!(λ::Vector{FT}, phot::Vector{FT}, E::Vector{FT}) where {FT<:AbstractFloat}

Compute and save the number of moles of photons, given
- `λ` Wave length in `[nm]`, converted to `[m]` by FAC
- `phot` Mole of photons (variable to save)
- `E` Joules of energy

---
# Examples
```julia
Es = zeros(2);
energy!.([500.0, 400.0], [0.1, 0.1], Es);
```
"""
energy!(λ::Vector{FT}, phot::Vector{FT}, E::Vector{FT}) where {FT<:AbstractFloat} = (
    E .= energy.(λ, phot);

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes made to this method
# General
#     2021-Jun-13: add method to save to provided 2rd variable
#
#######################################################################################################################################################################################################
"""

    energy!(λ::Vector{FT}, phot::Vector{FT}) where {FT<:AbstractFloat}

Compute and save the number of moles of photons, given
- `λ` Wave length in `[nm]`, converted to `[m]` by FAC
- `E` Joules of energy, will be converted to moles of photons

---
# Examples
```julia
phots = rand(2);
energy!.([500.0, 400.0], phots);
```
"""
energy!(λ::Vector{FT}, phot::Vector{FT}) where {FT<:AbstractFloat} = (
    phot ./= λ .* FT(FAC);

    return nothing
);
