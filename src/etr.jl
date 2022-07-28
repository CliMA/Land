#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: rename the function from leaf_ETR! to to photosystem_electron_transport! to be more specific
#     2022-Jan-18: use apar(ppar) and p_i as inputs rather than field from Leaf, and this allows for more modular operations
#     2022-Jan-18: use apar(ppar) rather than par as in Johnson and Berry's paper (they convert par to apar(ppar) later)
#     2022-Feb-07: add e_to_c calculation
#     2022-Feb-07: remove duplicated j (using j_pot is enough) for C4VJPModel
#     2022-Mar-01: save PSI J to psm._j_psi
#     2022-Mar-01: use η_c and η_l from psm (temperature corrected) rather than constant Η_C and Η_L
#     2022-Jul-01: add β to variable list to account for Vmax downregulation used in CLM5
#
#######################################################################################################################################################################################################
"""

    photosystem_electron_transport!(psm::C3CytochromeModel{FT}, rc::CytochromeReactionCenter{FT}, apar::FT, p_i::FT; β::FT = FT(1)) where {FT<:AbstractFloat}
    photosystem_electron_transport!(psm::C3VJPModel{FT}, rc::VJPReactionCenter{FT}, apar::FT, p_i::FT; β::FT = FT(1)) where {FT<:AbstractFloat}
    photosystem_electron_transport!(psm::C4VJPModel{FT}, rc::VJPReactionCenter{FT}, apar::FT, p_i::FT; β::FT = FT(1)) where {FT<:AbstractFloat}

Update the electron transport rates, given
- `psm` `C3CytochromeModel`, `C3VJPModel`, or `C4VJPModel` type C3 photosynthesis model
- `rc` `CytochromeReactionCenter` or `VJPReactionCenter` type photosynthesis system reaction center
- `apar` Absorbed photosynthetically active radiation in `μmol m⁻² s⁻¹`
- `p_i` Internal CO₂ partial pressure in `Pa`, used to compute e_to_c
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd

"""
function photosystem_electron_transport! end

photosystem_electron_transport!(psm::C3CytochromeModel{FT}, rc::CytochromeReactionCenter{FT}, apar::FT, p_i::FT; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack EFF_1, EFF_2 = psm;
    @unpack F_PSI, Φ_PSI_MAX = rc;

    psm._e_to_c = (p_i - psm._γ_star) / (EFF_1*p_i + EFF_2*psm._γ_star);
    psm._j_psi  = colimited_rate(β * psm._v_qmax, apar * F_PSI * Φ_PSI_MAX, psm.COLIMIT_J);
    psm._η      = 1 - psm._η_l / psm._η_c + (3*p_i + 7*psm._γ_star) / (EFF_1*p_i + EFF_2*psm._γ_star) / psm._η_c;
    psm._j_pot  = psm._j_psi / psm._η;

    return nothing
);

photosystem_electron_transport!(psm::C3VJPModel{FT}, rc::VJPReactionCenter{FT}, apar::FT, p_i::FT; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack EFF_1, EFF_2 = psm;
    @unpack F_PSII, Φ_PSII_MAX = rc;

    psm._e_to_c = (p_i - psm._γ_star) / (EFF_1*p_i + EFF_2*psm._γ_star);
    psm._j_pot  = F_PSII * Φ_PSII_MAX * apar;
    psm._j      = colimited_rate(psm._j_pot, β * psm._j_max, psm.COLIMIT_J);

    return nothing
);

photosystem_electron_transport!(psm::C4VJPModel{FT}, rc::VJPReactionCenter{FT}, apar::FT, p_i::FT; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    @unpack F_PSII, Φ_PSII_MAX = rc;

    psm._e_to_c = 1 / 6;
    psm._j_pot  = F_PSII * Φ_PSII_MAX * apar;

    return nothing
);
