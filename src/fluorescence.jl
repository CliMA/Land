#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: rename the function to photosystem_coefficients!
#     2022-Feb-07: remove the wrapper function method
#
#######################################################################################################################################################################################################
"""
This function updates the rate constants and coefficients in reaction center. Supported methods are

$(METHODLIST)

"""
function photosystem_coefficients! end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jan-14: unpack CONSTANT from the input variables only
#     2022-Jan-14: add function that operates PSM, PRC, and FLM directly so as to be more modular (reduce memory allocations)
#     2022-Jan-24: fix documentation
#     2022-Feb-07: use apar in fluorescence model (not used in this method)
#     2022-Feb-07: remove fluorescence model from input variables (in reaction center since ClimaCache v0.1.2)
# Bug fix
#     2022-Feb-24: a typo from "rc.ϕ_f  = rc.f_m′ / (1 - rc.ϕ_p);" to "rc.ϕ_f  = rc.f_m′ * (1 - rc.ϕ_p);"
#
#######################################################################################################################################################################################################
"""

    photosystem_coefficients!(psm::Union{C3VJPModel{FT}, C4VJPModel{FT}}, rc::VJPReactionCenter{FT}, apar::FT) where {FT<:AbstractFloat}

Update the rate constants and coefficients in reaction center, given
- `psm` `C3VJPModel` or `C4VJPModel` type photosynthesis model
- `rc` `VJPReactionCenter` type photosynthesis system reaction center
- `apar` Absorbed photosynthetically active radiation in `μmol m⁻² s⁻¹`
"""
photosystem_coefficients!(psm::Union{C3VJPModel{FT}, C4VJPModel{FT}}, rc::VJPReactionCenter{FT}, apar::FT) where {FT<:AbstractFloat} = (
    @unpack K_0, K_A, K_B = rc.FLM;
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
    rc.ϕ_f  = rc.f_m′ * (1 - rc.ϕ_p);

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
#     2022-Feb-07: remove fluorescence model from input variables
#     2022-Feb-07: use a_gross and j_pot rather than a series of j_p680 and j_p700
# Bug fix
#     2022-Feb-10: scale fluorescence quantum yield based on F_PSI and reabsorption factor
#     2022-Feb-10: _q1 needs to be multiply by η
# To do
#     TODO: add more calculations such as NPQ when the model is ready
#
#######################################################################################################################################################################################################
"""

    photosystem_coefficients!(psm::C3CytochromeModel{FT}, rc::CytochromeReactionCenter{FT}, apar::FT) where {FT<:AbstractFloat}

Update the rate constants and coefficients in reaction center, given
- `psm` `C3CytochromeModel` type photosynthesis model
- `rc` `CytochromeReactionCenter` type photosynthesis system reaction center
- `apar` Absorbed photosynthetically active radiation in `μmol m⁻² s⁻¹`
"""
photosystem_coefficients!(psm::C3CytochromeModel{FT}, rc::CytochromeReactionCenter{FT}, apar::FT) where {FT<:AbstractFloat} = (
    @unpack F_PSI, K_D, K_F, K_PSI, K_PSII, K_U, K_X, Φ_PSI_MAX = rc;

    # adapted from https://github.com/jenjohnson/johnson-berry-2021-pres/blob/main/scripts/model_fun.m
    _y2 = psm.a_gross / (psm.e_to_c * apar * (1 - F_PSI));
    _q2 = 1 - psm.j_pot * psm.η / psm.v_qmax;
    _q1 = psm.a_gross * psm.η / (psm.e_to_c * apar * F_PSI * Φ_PSI_MAX);

    # solve PSII K_N
    _k_n     = ( sqrt(K_PSII^2 * (_y2 - _q2)^2 + K_U^2 * _y2^2 + 2 * K_PSII * K_U * _y2 * (_y2 + _q2 - 2 * _y2 * _q2)) + K_PSII * (_q2 - _y2) + K_U * _y2 ) / (2 * _y2) - K_F - K_U - K_D;
    _k_sum_1 = K_D + K_F + K_U + _k_n;
    _k_sum_2 = K_D + K_F + K_U + _k_n + K_PSII;
    _k_sum_3 = K_D + K_F + K_PSI;
    _k_sum_4 = K_D + K_F + K_X;

    # compute PSII and PSI yeilds
    _ϕ_U2_a = _q2 * K_U    / _k_sum_2 + (1 - _q2) * K_U  / _k_sum_1;
    _ϕ_P2_a = _q2 * K_PSII / _k_sum_2 / (1 - _ϕ_U2_a);
    _ϕ_F2_a = (_q2 * K_F   / _k_sum_2 + (1 - _q2) * K_F  / _k_sum_1) / (1 - _ϕ_U2_a);
    _ϕ_F1_a = K_F / _k_sum_3 * _q1 + K_F / _k_sum_4 * (1 - _q1);

    # save the fluorescence and photosynthesis yields in reaction center
    rc.ϕ_f = _ϕ_F1_a * rc.ϵ_1 * F_PSI + _ϕ_F2_a * rc.ϵ_2 * (1 - F_PSI);
    rc.ϕ_p = _ϕ_P2_a * (1 - F_PSI);

    return nothing

    #=
    # some unused or unsaved variables
    #_ϕ_N2_a = (_q2 * _k_n  / _k_sum_2 + (1 - _q2) * _k_n / _k_sum_1) / (1 - _ϕ_U2_a);
    #_ϕ_D2_a = (_q2 * K_D   / _k_sum_2 + (1 - _q2) * K_D  / _k_sum_1) / (1 - _ϕ_U2_a);
    #_ϕ_P1_a = K_PSI / _k_sum_3 * _q1;
    #_ϕ_N1_a = K_X   / _k_sum_4 * (1 - _q1);
    #_ϕ_D1_a = K_D   / _k_sum_3 * _q1 + K_D / _k_sum_4 * (1 - _q1);

    # PAM measured fluorescence levels (Eqns. 38-42)
    #   N.B., hardcoding of a2(1) for dark-adapted value
    _tmp_1 = α_1 * K_F * ϵ_1;
    _tmp_2 = α_2 * K_F * ϵ_2;
    _Fm_a  = _tmp_1 / _k_sum_4 + _tmp_2 / (K_D + K_F);
    _Fo_a  = _tmp_1 / _k_sum_3 + _tmp_2 / (K_D + K_F + K_PSII);
    _Fmp_a = _tmp_1 / _k_sum_4 + _tmp_2 / (K_D + K_F + _k_n);
    _Fop_a = _tmp_1 / _k_sum_3 + _tmp_2 / (K_D + K_F + K_PSII + _k_n);
    _Fs_a  = α_1 * _ϕ_F1_a * ϵ_1 + α_2 * _ϕ_F2_a * ϵ_2;

    # PAM indices used in plotter_forward_fun.m
    _PAM1_a = 1 - _Fs_a / _Fmp_a; # ϕ_P
    _PAM2_a = _Fs_a * (1 / _Fmp_a - 1/_Fm_a); # ϕ_N
    _PAM3_a = _Fs_a / _Fm_a; # ϕ_D + ϕ_F
    =#
);
