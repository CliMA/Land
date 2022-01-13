###############################################################################
#
# Update the fluorescence in the leaf
#
###############################################################################
"""
    leaf_fluorescence!(
                fluo_set::CytoFluoParaSet{FT},
                leaf::Leaf{FT}
    ) where {FT<:AbstractFloat}
    leaf_fluorescence!(
                fluo_set::FluoParaSet{FT},
                leaf::Leaf{FT}
    ) where {FT<:AbstractFloat}

Compute fluorescence yield, Kr, Ks, and Kp for leaf, given
- `fluo_set` [`FluoParaSet`](@ref) type parameter set
- `leaf` [`Leaf`](@ref) struct
"""
function leaf_fluorescence!(
            fluo_set::CytoFluoParaSet{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    @unpack APAR, C_b₆f, J_P680_a, J_P700_a, J_P700_j, K_D1, K_F1, K_N1, K_P1,
            K_P2, K_U2, k_q, φ_P1_max, α_1, α_2, ϵ_1, ϵ_2 = leaf;

    # adapted from https://github.com/jenjohnson/johnson-berry-2021-pres/
    #                      scripts/model_fun.m
    # TODO PAR or APAR?
    # primary fluorescence parameters
    _b₆f_a  = J_P700_j / k_q;
    _φ_P1_a = J_P700_a / (APAR * α_1);
    _q1_a   = _φ_P1_a  / φ_P1_max;
    _φ_P2_a = J_P680_a / (APAR * α_2);
    _q2_a   = 1 - _b₆f_a / C_b₆f;

    # N.B., rearrange Eqn. 25a to solve for _kn_2_a
    _kn_2_a = ( sqrt(K_P2^2 * _φ_P2_a^2 -
                     2 * K_P2^2 * _φ_P2_a * _q2_a +
                     K_P2^2 * _q2_a^2 -
                     4 * K_P2 * K_U2 * _φ_P2_a^2 * _q2_a +
                     2 * K_P2 * K_U2 * _φ_P2_a^2 +
                     2 * K_P2 * K_U2 * _φ_P2_a * _q2_a +
                     K_U2^2 * _φ_P2_a^2) -
                K_P2 * _φ_P2_a +
                K_U2 * _φ_P2_a +
                K_P2 * _q2_a
              ) / (2 * _φ_P2_a) - K_F1 - K_U2 - K_D1;

    # Photosystem II (Eqns. 23a-23e and 25a-25d)
    _k_sum_1 = _kn_2_a + K_D1 + K_F1 + K_U2;
    _k_sum_2 = _kn_2_a + K_D1 + K_F1 + K_U2 + K_P2;
    _φ_P2_a  = _q2_a * K_P2    / _k_sum_2;
    _φ_N2_a  = _q2_a * _kn_2_a / _k_sum_2 + (1 - _q2_a) * _kn_2_a / _k_sum_1;
    _φ_D2_a  = _q2_a * K_D1    / _k_sum_2 + (1 - _q2_a) * K_D1    / _k_sum_1;
    _φ_F2_a  = _q2_a * K_F1    / _k_sum_2 + (1 - _q2_a) * K_F1    / _k_sum_1;
    _φ_U2_a  = _q2_a * K_U2    / _k_sum_2 + (1 - _q2_a) * K_U2    / _k_sum_1;
    _φ_P2_a /= 1 - _φ_U2_a;
    _φ_N2_a /= 1 - _φ_U2_a;
    _φ_D2_a /= 1 - _φ_U2_a;
    _φ_F2_a /= 1 - _φ_U2_a;
    leaf.φ = _φ_F2_a;

    #=
    # For Photosystem I (Eqns. 19a-19d)
    _k_sum_3 = K_P1 + K_D1 + K_F1;
    _k_sum_4 = K_N1 + K_D1 + K_F1;
    _φ_P1_a  = _q1_a * K_P1 / _k_sum_3;
    _φ_N1_a  = (1 - _q1_a) * K_N1 / _k_sum_4;
    _φ_D1_a  = _q1_a * K_D1 / _k_sum_3 + (1 - _q1_a) * K_D1 / _k_sum_4;
    _φ_F1_a  = _q1_a * K_F1 / _k_sum_3 + (1 - _q1_a) * K_F1 / _k_sum_4;

    # PAM measured fluorescence levels (Eqns. 38-42)
    #   N.B., hardcoding of a2(1) for dark-adapted value
    _tmp_1 = α_1 * K_F1 * ϵ_1;
    _tmp_2 = α_2 * K_F1 * ϵ_2;
    _Fm_a  = _tmp_1 / _k_sum_4 + _tmp_2 / (K_D1 + K_F1);
    _Fo_a  = _tmp_1 / _k_sum_3 + _tmp_2 / (K_D1 + K_F1 + K_P2);
    _Fmp_a = _tmp_1 / _k_sum_4 + _tmp_2 / (K_D1 + K_F1 + _kn_2_a);
    _Fop_a = _tmp_1 / _k_sum_3 + _tmp_2 / (K_D1 + K_F1 + K_P2 + _kn_2_a);
    _Fs_a  = α_1 * _φ_F1_a * ϵ_1 + α_2 * _φ_F2_a * ϵ_2;

    # PAM indices used in plotter_forward_fun.m
    _PAM1_a = 1 - _Fs_a / _Fmp_a; # φ_P
    _PAM2_a = _Fs_a * (1 / _Fmp_a - 1/_Fm_a); # φ_N
    _PAM3_a = _Fs_a / _Fm_a; # φ_D + φ_F
    =#

    return nothing
end




function leaf_fluorescence!(
            fluo_set::FluoParaSet{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    @unpack Ag, Kd, Kf, Kp_max, maxPSII = leaf;
    @unpack Kr1, Kr2, Kr3 = fluo_set;

    # Actual effective ETR:
    leaf.Ja  = max(0, Ag / leaf.e2c);

    # Effective photochemical yield:
    if leaf.Ja <= 0
        _φ   = maxPSII;
    else
        _φ   = maxPSII*leaf.Ja/leaf.J_pot;
    end

    leaf.φ   = min(1/maxPSII, _φ);
    # degree of light saturation: 'x' (van der Tol e.Ap. 2014)
    x        = max(0,  1-leaf.φ/maxPSII);

    # Max PSII rate constant, x_α = exp(log(x)*Kr2);
    x_α      = x ^ Kr2;
    leaf.Kr  = Kr1 * (1 + Kr3)* x_α / (Kr3 + x_α);
    leaf.Kp  = max(0, leaf.φ*(Kf+Kd+leaf.Kr)/(1-leaf.φ));

    leaf.Fo  = Kf / (Kf+Kp_max+Kd                );
    leaf.Fo′ = Kf / (Kf+Kp_max+Kd+leaf.Kr+leaf.Ks);
    leaf.Fm  = Kf / (Kf       +Kd                );
    leaf.Fm′ = Kf / (Kf       +Kd+leaf.Kr+leaf.Ks);
    leaf.φs  = leaf.Fm′ * (1 - leaf.φ);

    #
    #
    # TODO should these change with Ks?
    #
    #
    # leaf.eta  = leaf.ϕs/leaf.Fo
    # don't need this anymore
    # better to use ϕs directly for SIF as Fo is not always fqe=0.01
    leaf.qQ  = 1 - (leaf.φs-leaf.Fo′) / (leaf.Fm-leaf.Fo′);
    leaf.qE  = 1 - (leaf.Fm-leaf.Fo′) / (leaf.Fm′-leaf.Fo);
    leaf.NPQ = leaf.Kr / (Kf + Kd);

    return nothing
end
