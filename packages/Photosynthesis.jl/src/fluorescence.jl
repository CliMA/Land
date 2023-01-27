#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: rename the function to photosystem_coefficients!
#     2022-Jan-14: add function that operates PSM, PRC, and FLM directly so as to be more modular (reduce memory allocations)
#     2022-Feb-07: remove the wrapper function method
#     2022-Feb-07: use ppar in fluorescence model (not used in this method)
#     2022-Feb-07: remove fluorescence model from input variables (in reaction center since ClimaCache v0.1.2)
#     2022-Feb-07: add support for Johnson and Berry (2021) model
#     2022-Feb-07: use a_gross and j_pot rather than a series of j_p680 and j_p700
#     2022-Feb-10: scale fluorescence quantum yield based on F_PSI and reabsorption factor
#     2022-Feb-10: _q1 needs to be multiply by η
#     2022-Mar-04: add support to sustained NPQ
#     2022-Mar-04: use the weighted yield for photosynthesis
#     2022-Jul-01: add β to variable list to account for Vmax downregulation used in CLM5
# Bug fix
#     2022-Feb-24: a typo from "rc.ϕ_f  = rc.f_m′ / (1 - rc.ϕ_p);" to "rc.ϕ_f  = rc.f_m′ * (1 - rc.ϕ_p);"
#     2022-Feb-28: psm.e_to_c is recalculated based on analytically resolving leaf.p_CO₂_i from leaf.g_CO₂, this psm.e_to_c used to be calculated as psm.a_j / psm.j (a_j here is not p_CO₂_i based)
#                  note here that in CliMA v0.1, this e_to_c is not updated properly.
# To do
#     TODO: add more calculations such as NPQ when the model is ready
#
#######################################################################################################################################################################################################
"""

    photosystem_coefficients!(psm::C3CytochromeModel{FT}, rc::CytochromeReactionCenter{FT}, ppar::FT; β::FT = FT(1)) where {FT<:AbstractFloat}
    photosystem_coefficients!(psm::Union{C3VJPModel{FT}, C4VJPModel{FT}}, rc::VJPReactionCenter{FT}, ppar::FT; β::FT = FT(1)) where {FT<:AbstractFloat}

Update the rate constants and coefficients in reaction center, given
- `psm` `C3CytochromeModel`, `C3VJPModel`, or `C4VJPModel` type photosynthesis model
- `rc` `CytochromeReactionCenter` or `VJPReactionCenter` type photosynthesis system reaction center
- `ppar` Absorbed photosynthetically active radiation in `μmol m⁻² s⁻¹`
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd

"""
function photosystem_coefficients! end

photosystem_coefficients!(psm::C3CytochromeModel{FT}, rc::CytochromeReactionCenter{FT}, ppar::FT; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    (; F_PSI, K_D, K_F, K_PSI, K_PSII, K_U, K_X, Φ_PSI_MAX) = rc;

    # adapted from https://github.com/jenjohnson/johnson-berry-2021-pres/blob/main/scripts/model_fun.m
    _ϕ_P1_a = psm.a_gross * psm._η / (psm._e_to_c * ppar * F_PSI);
    _ϕ_P2_a = psm.a_gross / (psm._e_to_c * ppar * (1 - F_PSI));
    _q1     = _ϕ_P1_a / Φ_PSI_MAX;
    _q2     = 1 - psm._j_psi / (β * psm._v_qmax);

    # solve PSII K_N
    _k_sum_na = _ϕ_P2_a;
    _k_sum_nb = -1 * (K_U * _ϕ_P2_a + K_PSII * (_q2 - _ϕ_P2_a));
    _k_sum_nc = -1 * (_ϕ_P2_a * (1 - _q2) * K_U * K_PSII);
    _k_sum    = upper_quadratic(_k_sum_na, _k_sum_nb, _k_sum_nc);
    _k_n      = _k_sum - K_F - K_U - K_D;

    # compute PSII and PSI yeilds
    _k_sum_1 = K_D + K_F + K_U + _k_n;
    _k_sum_2 = K_D + K_F + K_U + _k_n + K_PSII;
    _k_sum_3 = K_D + K_F + K_PSI;
    _k_sum_4 = K_D + K_F + K_X;
    _ϕ_U2_a  =  _q2 * K_U / _k_sum_2 + (1 - _q2) * K_U  / _k_sum_1;
    _ϕ_F2_a  = (_q2 * K_F / _k_sum_2 + (1 - _q2) * K_F  / _k_sum_1) / (1 - _ϕ_U2_a);
    _ϕ_F1_a  = K_F / _k_sum_3 * _q1 + K_F / _k_sum_4 * (1 - _q1);

    # save the weighted fluorescence and photosynthesis yields in reaction center
    rc.ϕ_f = _ϕ_F1_a * rc.ϵ_1 * F_PSI + _ϕ_F2_a * rc.ϵ_2 * (1 - F_PSI);
    rc.ϕ_p = _ϕ_P1_a * F_PSI + _ϕ_P2_a * (1 - F_PSI);

    return nothing

    #=
    # some unused or unsaved variables
    _ϕ_N2_a = (_q2 * _k_n  / _k_sum_2 + (1 - _q2) * _k_n / _k_sum_1) / (1 - _ϕ_U2_a);
    _ϕ_D2_a = (_q2 * K_D   / _k_sum_2 + (1 - _q2) * K_D  / _k_sum_1) / (1 - _ϕ_U2_a);
    _ϕ_P1_a = K_PSI / _k_sum_3 * _q1;
    _ϕ_N1_a = K_X   / _k_sum_4 * (1 - _q1);
    _ϕ_D1_a = K_D   / _k_sum_3 * _q1 + K_D / _k_sum_4 * (1 - _q1);

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

photosystem_coefficients!(psm::Union{C3VJPModel{FT}, C4VJPModel{FT}}, rc::VJPReactionCenter{FT}, ppar::FT; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    (; K_0, K_A, K_B) = rc.FLM;
    (; F_PSII, K_D, K_F, K_P_MAX, Φ_PSII_MAX) = rc;

    # calculate photochemical yield
    rc.ϕ_p = psm.a_gross / (psm._e_to_c * F_PSII * ppar);

    # calculate rate constants
    _x            = max(0, 1 - rc.ϕ_p / Φ_PSII_MAX);
    _xᵅ           = _x ^ K_A;
    rc._k_npq_rev = K_0 * (1 + K_B) * _xᵅ / (K_B + _xᵅ);
    rc._k_p       = max(0, rc.ϕ_p * (K_F + K_D + rc._k_npq_rev + rc.k_npq_sus) / (1 - rc.ϕ_p) );

    # TODO: whether to consider sustained K_N in the calculations of f_o and f_m
    # rc._f_o  = K_F / (K_F + K_P_MAX + K_D + rc.k_npq_sus);
    # rc._f_o′ = K_F / (K_F + K_P_MAX + K_D + rc.k_npq_sus + rc._k_npq_rev);
    # rc._f_m  = K_F / (K_F + K_D + rc.k_npq_sus);
    # rc._f_m′ = K_F / (K_F + K_D + rc.k_npq_sus + rc._k_npq_rev);

    # calculate fluorescence quantum yield
    rc._f_o  = K_F / (K_F + K_P_MAX + K_D);
    rc._f_o′ = K_F / (K_F + K_P_MAX + K_D + rc._k_npq_rev + rc.k_npq_sus);
    rc._f_m  = K_F / (K_F + K_D);
    rc._f_m′ = K_F / (K_F + K_D + rc._k_npq_rev + rc.k_npq_sus);
    rc.ϕ_f   = rc._f_m′ * (1 - rc.ϕ_p);

    # TODO: if K_N is used above, do we need to recalculate _npq
    # rc._npq = (rc._k_npq_rev + rc.k_npq_sus) / (K_F + K_D + rc.k_npq_sus);

    # calculate quenching rates
    rc._q_e = 1 - (rc._f_m - rc._f_o′) / (rc._f_m′ - rc._f_o);
    rc._q_p = 1 - (rc.ϕ_f - rc._f_o′) / (rc._f_m - rc._f_o′);
    rc._npq = (rc._k_npq_rev + rc.k_npq_sus) / (K_F + K_D);

    return nothing
);
