###############################################################################
#
# Calculate the electron transport rate
#
###############################################################################
"""
    leaf_ETR!(photo_set::AbstractPhotoModelParaSet{FT},
              leaf::Leaf{FT}
    ) where {FT<:AbstractFloat}

Update the electron transport variables in the leaf struct, given
- `photo_set` [`C3ParaSet`](@ref) or [`C4ParaSet`](@ref) type struct
- `leaf` [`Leaf`](@ref) type struct
"""
function leaf_ETR!(
            photo_set::C3Cytochrome{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    @unpack APAR, K_P1, Vqmax, n_C, n_L, p_i, α_1, Γ_star, φ_P1_max = leaf;
    @unpack Eff_1, Eff_2 = photo_set;

    # TODO note to run this function before rubisco_limited_rate!
    # TODO Is Vqmax necessary here, or we use Jmax as well
    # TODO Johnson and Berry used PAR here
    leaf.J_P700_j = Vqmax * APAR / (Vqmax / α_1 / φ_P1_max + APAR);
    leaf.η = 1 - n_L / n_C + (3*p_i + 7*Γ_star) / (Eff_1*p_i + Eff_2*Γ_star) / n_C;
    leaf.J_P680_j = leaf.J_P700_j / leaf.η;

    leaf.J_pot = leaf.J_P680_j;
    leaf.J = leaf.J_P680_j;

    return nothing
end




function leaf_ETR!(
            photo_set::C3ParaSet{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    @unpack APAR, maxPSII, Jmax, PSII_frac = leaf;

    _Jp = PSII_frac * maxPSII * APAR;
    _J  = min(_Jp, Jmax);

    leaf.J_pot = _Jp;
    leaf.J     = _J;

    return nothing
end




function leaf_ETR!(
            photo_set::C4ParaSet{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    @unpack APAR, maxPSII, PSII_frac = leaf

    _Jp = PSII_frac * maxPSII * APAR;

    leaf.J_pot = _Jp;
    leaf.J     = _Jp;

    return nothing
end
