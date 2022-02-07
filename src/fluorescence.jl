#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: rename the function to photosystem_coefficients!
#
#######################################################################################################################################################################################################
"""
This function updates the rate constants and coefficients in reaction center. Supported methods are

$(METHODLIST)

"""
function photosystem_coefficients! end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: add function that for simple function call directly on Leaf
#     2022-Jan-24: fix documentation
#
#######################################################################################################################################################################################################
"""

    photosystem_coefficients!(leaf::Leaf{FT}) where {FT<:AbstractFloat}

Update the rate constants and coefficients in reaction center, given
- `leaf` `Leaf` type structure that stores biophysical, reaction center, and photosynthesis model structures
"""
photosystem_coefficients!(leaf::Leaf{FT}) where {FT<:AbstractFloat} = (
    photosystem_coefficients!(leaf.PSM, leaf.PRC, leaf.FLM);

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jan-14: unpack CONSTANT from the input variables only
#     2022-Jan-14: add function that operates PSM, PRC, and FLM directly so as to be more modular (reduce memory allocations)
#     2022-Jan-24: fix documentation
#
#######################################################################################################################################################################################################
"""

    photosystem_coefficients!(psm::Union{C3VJPModel{FT}, C4VJPModel{FT}}, rc::VJPReactionCenter{FT}, vdt::VanDerTolFluorescenceModel{FT}) where {FT<:AbstractFloat}

Update the rate constants and coefficients in reaction center, given
- `psm` `C3VJPModel` or `C4VJPModel` type photosynthesis model
- `rc` `VJPReactionCenter` type photosynthesis system reaction center
- `vdt` van der Tol et al. (2013) fluorescence model
"""
photosystem_coefficients!(psm::Union{C3VJPModel{FT}, C4VJPModel{FT}}, rc::VJPReactionCenter{FT}, vdt::VanDerTolFluorescenceModel{FT}) where {FT<:AbstractFloat} = (
    @unpack K_0, K_A, K_B = vdt;
    @unpack K_D, K_F, K_P_MAX, Φ_PSII_MAX = rc;

    # calculate photochemical yield
    rc.ϕ_p = psm.a_gross / psm.e_to_c / psm.j_pot * Φ_PSII_MAX;

    # calculate rate constants
    _x           = max(0, 1 - rc.ϕ_p / Φ_PSII_MAX);
    _xᵅ          = _x ^ K_A;
    rc.k_npq_rev = K_0 * (1 + K_B) * _xᵅ / (K_B + _xᵅ);
    rc.k_p       = max(0, rc.ϕ_p * (K_F + K_D + rc.k_npq_rev) / (1 - rc.ϕ_p) );

    # calculate fluorescence quantum yield
    rc.f_o  = K_F / (K_F + K_P_MAX + K_D);
    rc.f_o′ = K_F / (K_F + K_P_MAX + K_D + rc.k_npq_rev + rc.k_npq_sus);
    rc.f_m  = K_F / (K_F + K_D);
    rc.f_m′ = K_F / (K_F + K_D + rc.k_npq_rev + rc.k_npq_sus);
    rc.ϕ_f  = rc.f_m′ / (1 - rc.ϕ_p);

    # calculate quenching rates
    rc.q_e = 1 - (rc.f_m - rc.f_o′) / (rc.f_m′ - rc.f_o);
    rc.q_p = 1 - (rc.ϕ_f - rc.f_o′) / (rc.f_m - rc.f_o′);
    rc.npq = rc.k_npq_rev / (K_F + K_D);

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Feb-07: add support for Johnson and Berry (2021) model
#
#######################################################################################################################################################################################################
"""

    photosystem_coefficients!(psm::C3CytochromeModel{FT}, rc::CytochromeReactionCenter{FT}, vdt::CytochromeFluorescenceModel{FT}) where {FT<:AbstractFloat}

Update the rate constants and coefficients in reaction center, given
- `psm` `C3CytochromeModel` type photosynthesis model
- `rc` `CytochromeReactionCenter` type photosynthesis system reaction center
- `cfm` Johnson and Berry (2021) fluorescence model
"""
photosystem_coefficients!(psm::C3CytochromeModel{FT}, rc::CytochromeReactionCenter{FT}, cfm::CytochromeFluorescenceModel{FT}) where {FT<:AbstractFloat} = (
    return nothing
)
