#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: rename the function from leaf_ETR! to to photosystem_electron_transport! to be more specific
#
#######################################################################################################################################################################################################
"""
This function updates the electron transport rates in photosynthesis reaction center, and thus to calculate photosynthetic rate. Supported methods are

$(METHODLIST)

"""
function photosystem_electron_transport! end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: use apar as an input rather than field from Leaf, and this allows for more modular operations
#     2022-Jan-14: remove examples from the doc, as this function is not meant to be public
#     2022-Jan-14: unpack CONSTANT from the input variables only
#     2022-Jan-18: add input variable p_i for modularity
#     2022-Jan-24: fix documentation
#     2022-Feb-28: move e_to_c calculation back into this function to get aligned with C3CytochromeModel at GCO₂Mode analytically
#
#######################################################################################################################################################################################################
"""

    photosystem_electron_transport!(psm::C3VJPModel{FT}, rc::VJPReactionCenter{FT}, apar::FT, p_i::FT) where {FT<:AbstractFloat}

Update the electron transport rates, given
- `psm` `C3VJPModel` type C3 photosynthesis model
- `rc` `VJPReactionCenter` type photosynthesis system reaction center
- `apar` Absorbed photosynthetically active radiation in `μmol m⁻² s⁻¹`
- `p_i` Internal CO₂ partial pressure in `Pa`, used to compute e_to_c
"""
photosystem_electron_transport!(psm::C3VJPModel{FT}, rc::VJPReactionCenter{FT}, apar::FT, p_i::FT) where {FT<:AbstractFloat} = (
    @unpack EFF_1, EFF_2 = psm;
    @unpack F_PSII, Φ_PSII_MAX = rc;

    psm.e_to_c = (p_i - psm.γ_star) / (EFF_1*p_i + EFF_2*psm.γ_star);
    psm.j_pot  = F_PSII * Φ_PSII_MAX * apar;
    psm.j      = colimited_rate(psm.j_pot, psm.j_max, psm.COLIMIT_J);

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-18: use apar and p_i as inputs rather than field from Leaf, and this allows for more modular operations
#     2022-Jan-18: unpack CONSTANT from the input variables only
#     2022-Jan-18: use apar rather than par as in Johnson and Berry's paper
#     2022-Jan-24: fix documentation
#     2022-Feb-07: use correct field names in C3CytochromeModel
#     2022-Feb-07: remove j_pot and j definitions
#     2022-Feb-07: add e_to_c calculation
#     2022-Feb-07: add j_pot calculation back to work with ClimaCache v0.1.2
#     2022-Feb-07: move e_to_c calculation back to light_limited_rate!
#     2022-Feb-28: move e_to_c calculation back into this function to get aligned with C3CytochromeModel at GCO₂Mode analytically
#     2022-Mar-01: save PSI J to psm.j_psi
#     2022-Mar-01: use η_c and η_l from psm (temperature corrected) rather than constant Η_C and Η_L
#
#######################################################################################################################################################################################################
"""

    photosystem_electron_transport!(psm::C3CytochromeModel{FT}, rc::CytochromeReactionCenter{FT}, apar::FT, p_i::FT) where {FT<:AbstractFloat}

Update the electron transport rates, given
- `psm` `C3CytochromeModel` type C4 photosynthesis model
- `rc` `CytochromeReactionCenter` type photosynthesis system reaction center
- `apar` Absorbed photosynthetically active radiation in `μmol m⁻² s⁻¹`
- `p_i` Internal CO₂ partial pressure in `Pa`
"""
photosystem_electron_transport!(psm::C3CytochromeModel{FT}, rc::CytochromeReactionCenter{FT}, apar::FT, p_i::FT) where {FT<:AbstractFloat} = (
    @unpack EFF_1, EFF_2 = psm;
    @unpack F_PSI, Φ_PSI_MAX = rc;

    psm.e_to_c = (p_i - psm.γ_star) / (EFF_1*p_i + EFF_2*psm.γ_star);
    psm.j_psi  = colimited_rate(psm.v_qmax, apar * F_PSI * Φ_PSI_MAX, psm.COLIMIT_J);
    psm.η      = 1 - psm.η_l / psm.η_c + (3*p_i + 7*psm.γ_star) / (EFF_1*p_i + EFF_2*psm.γ_star) / psm.η_c;
    psm.j_pot  = psm.j_psi / psm.η;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: use apar as an input rather than field from Leaf, and this allows for more modular operations
#     2022-Jan-14: remove examples from the doc, as this function is not meant to be public
#     2022-Jan-14: unpack CONSTANT from the input variables only
#     2022-Jan-18: add input variable p_i for modularity
#     2022-Jan-24: fix documentation
#     2022-Feb-07: add e_to_c definition
#     2022-Feb-07: remove duplicated j (using j_pot is enough)
#     2022-Feb-28: move e_to_c calculation back into this function to get aligned with C3CytochromeModel at GCO₂Mode analytically
#
#######################################################################################################################################################################################################
"""

    photosystem_electron_transport!(psm::C4VJPModel{FT}, rc::VJPReactionCenter{FT}, apar::FT, p_i::FT) where {FT<:AbstractFloat}

Update the electron transport rates, given
- `psm` `C4VJPModel` type C4 photosynthesis model
- `rc` `VJPReactionCenter` type photosynthesis system reaction center
- `apar` Absorbed photosynthetically active radiation in `μmol m⁻² s⁻¹`
- `p_i` Internal CO₂ partial pressure in `Pa`, not used in this method
"""
photosystem_electron_transport!(psm::C4VJPModel{FT}, rc::VJPReactionCenter{FT}, apar::FT, p_i::FT) where {FT<:AbstractFloat} = (
    @unpack F_PSII, Φ_PSII_MAX = rc;

    psm.e_to_c = 1 / 6;
    psm.j_pot  = F_PSII * Φ_PSII_MAX * apar;

    return nothing
);
