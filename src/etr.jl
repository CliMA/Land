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
#
#######################################################################################################################################################################################################
"""

    photosystem_electron_transport!(psm::C3VJPModel{FT}, rc::VJPReactionCenter{FT}, apar::FT, p_i::FT = FT(0)) where {FT<:AbstractFloat}

Update the electron transport rates, given
- `psm` `C3VJPModel` type C3 photosynthesis model
- `rc` `VJPReactionCenter` type photosynthesis system reaction center
- `apar` Absorbed photosynthetically active radiation in `μmol m⁻² s⁻¹`
- `p_i` Internal CO₂ partial pressure in `Pa`, not used in this method
"""
photosystem_electron_transport!(psm::C3VJPModel{FT}, rc::VJPReactionCenter{FT}, apar::FT, p_i::FT = FT(0)) where {FT<:AbstractFloat} = (
    @unpack F_PSII, Φ_PSII_MAX = rc;

    psm.j_pot = F_PSII * Φ_PSII_MAX * apar;
    psm.j     = min(psm.j_pot, psm.j_max);

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
#
#######################################################################################################################################################################################################
"""

    photosystem_electron_transport!(psm::C3CytochromeModel{FT}, rc::CytochromeReactionCenter{FT}, apar::FT, p_i::FT) where {FT<:AbstractFloat}

Update the electron transport rates, given
- `psm` `C3CytochromeModel` type C4 photosynthesis model
- `rc` `CytochromeReactionCenter` type photosynthesis system reaction center
- `apar` Absorbed photosynthetically active radiation in `μmol m⁻² s⁻¹`
- `p_i` Internal CO₂ partial pressure in `Pa`, not used in this method
"""
photosystem_electron_transport!(psm::C3CytochromeModel{FT}, rc::CytochromeReactionCenter{FT}, apar::FT, p_i::FT) where {FT<:AbstractFloat} = (
    @unpack EFF_1, EFF_2 = psm;
    @unpack F_PSI, Η_C, Η_L, Φ_PSI_MAX = rc;

    psm.j_p700_j = psm.v_qmax * apar * F_PSI * Φ_PSI_MAX / (psm.v_qmax + apar * F_PSI * Φ_PSI_MAX);
    psm.η        = 1 - Η_L / Η_C + (3*p_i + 7*psm.γ_star) / (EFF_1*p_i + EFF_2*psm.γ_star) / Η_C;
    psm.j_p680_j = psm.j_p700_j / psm.η;
    psm.j_pot    = psm.j_p680_j;
    psm.j        = psm.j_p680_j;

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
#
#######################################################################################################################################################################################################
"""

    photosystem_electron_transport!(psm::C4VJPModel{FT}, rc::VJPReactionCenter{FT}, apar::FT, p_i::FT = FT(0)) where {FT<:AbstractFloat}

Update the electron transport rates, given
- `psm` `C4VJPModel` type C4 photosynthesis model
- `rc` `VJPReactionCenter` type photosynthesis system reaction center
- `apar` Absorbed photosynthetically active radiation in `μmol m⁻² s⁻¹`
- `p_i` Internal CO₂ partial pressure in `Pa`, not used in this method
"""
photosystem_electron_transport!(psm::C4VJPModel{FT}, rc::VJPReactionCenter{FT}, apar::FT, p_i::FT = FT(0)) where {FT<:AbstractFloat} = (
    @unpack F_PSII, Φ_PSII_MAX = rc;

    psm.j_pot = F_PSII * Φ_PSII_MAX * apar;
    psm.j     = psm.j_pot;

    return nothing
);
