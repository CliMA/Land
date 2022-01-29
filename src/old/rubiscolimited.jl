#=
function rubisco_limited_rate!(
            photo_set::C3Cytochrome{FT},
            leaf::Leaf{FT}
) where {FT<:AbstractFloat}
    @unpack Km, p_i, Vcmax, Γ_star = leaf;
    @unpack Eff_1, Eff_2 = photo_set;

    leaf.Ac = Vcmax * (p_i - Γ_star) / (p_i + Km);
    leaf.J_P680_c = leaf.Ac * (Eff_1*p_i + Eff_2*Γ_star) / (p_i - Γ_star);
    leaf.J_P700_c = leaf.J_P680_c * leaf.η;

    return nothing
end
=#
