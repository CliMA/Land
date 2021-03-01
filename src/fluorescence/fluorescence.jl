###############################################################################
#
# Update the fluorescence in the leaf
#
###############################################################################
"""
    leaf_fluorescence!(
                fluo_set::FluoParaSet{FT},
                leaf::Leaf{FT}
    ) where {FT<:AbstractFloat}

Compute fluorescence yield, Kr, Ks, and Kp for leaf, given
- `fluo_set` [`FluoParaSet`](@ref) type parameter set
- `leaf` [`Leaf`](@ref) struct
"""
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
    leaf.ϕs  = leaf.Fm′ * (1 - leaf.φ);

    #
    #
    # TODO should these change with Ks?
    #
    #
    # leaf.eta  = leaf.ϕs/leaf.Fo
    # don't need this anymore
    # better to use ϕs directly for SIF as Fo is not always fqe=0.01
    leaf.qQ  = 1 - (leaf.ϕs-leaf.Fo′) / (leaf.Fm-leaf.Fo′);
    leaf.qE  = 1 - (leaf.Fm-leaf.Fo′) / (leaf.Fm′-leaf.Fo);
    leaf.NPQ = leaf.Kr / (Kf + Kd);

    return nothing
end
