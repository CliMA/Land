"""
This function updates the rate constants and coefficients in reaction center. Supported methods are

$(METHODLIST)
"""
function leaf_fluorescence! end


"""
    leaf_fluorescence!(ps::Union{C3VJPModel{FT}, C4VJPModel{FT}}, rc::PhotosynthesisReactionCenter{FT}, vdt::VanDerTolFluorescenceModel{FT}) where {FT<:AbstractFloat}

Update the rate constants and coefficients in reaction center, given
- `ps` `C3VJPModel` or `C4VJPModel` type photosynthesis model
- `rc` `PhotosynthesisReactionCenter` reaction center for rate constants and coefficients
- `vdt` van der Tol et al. (2013) fluorescence model
"""
leaf_fluorescence!(ps::Union{C3VJPModel{FT}, C4VJPModel{FT}}, rc::PhotosynthesisReactionCenter{FT}, vdt::VanDerTolFluorescenceModel{FT}) where {FT<:AbstractFloat} = (
    @unpack K_0, K_A, K_B = vdt;
    @unpack K_D, K_F, K_P_MAX = rc;

    # calculate photochemical yield
    rc.ϕ_p = ps.a_gross / ps.e_to_c / ps.j_pot * rc.Φ_PSII_MAX;

    # calculate rate constants
    _x           = max(0, 1 - rc.ϕ_p / rc.Φ_PSII_MAX);
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


"""
    leaf_fluorescence!(leaf::Leaf{FT}) where {FT<:AbstractFloat}

Update the rate constants and coefficients in reaction center, given
- `leaf` `Leaf` type structure that stores biophysical, reaction center, and photosynthesis model structures
"""
leaf_fluorescence!(leaf::Leaf{FT}) where {FT<:AbstractFloat} = (
    leaf_fluorescence!(leaf.PSM, leaf.PRC, leaf.FLM);

    return nothing
);
