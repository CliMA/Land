###############################################################################
#
# Update the fluorescence in the leaf
#
###############################################################################
"""
    leaf_fluorescence!(fluo_set::FluoParaSet{FT}, leaf::Leaf{FT})

Compute fluorescence yield, Kn, and Kp for leaf, given
- `fluo_set` [`FluoParaSet`](@ref) type parameter set
- `leaf` [`Leaf`](@ref) struct
"""
function leaf_fluorescence!(
            fluo_set::FluoParaSet{FT},
            leaf::Leaf{FT}
            ) where {FT<:AbstractFloat}
    @unpack Ag, Kd, Kf, maxPSII = leaf;
    @unpack Kn1, Kn2, Kn3 = fluo_set;

    # Actual effective ETR:
    leaf.Ja  = max(0, Ag / leaf.CO₂_per_electron);

    # Effective photochemical yield:
    if leaf.Ja <= 0
        _φ   = maxPSII;
    else
        _φ   = maxPSII*leaf.Ja/leaf.J_pot;
    end

    leaf.φ   = min(1/maxPSII, _φ);
    # degree of light saturation: 'x' (van der Tol e.Ap. 2014)
    x        = max(0,  1-leaf.φ/maxPSII);

    # Max PSII rate constant
    Kp_max   = FT(4.0);

    x_alpha  = exp(log(x)*Kn2);
    #println(x_alpha)

    leaf.Kn  = Kn1 * (1+Kn3)* x_alpha/(Kn3 + x_alpha);
    leaf.Kp  = max(0,leaf.φ*(Kf+Kd+leaf.Kn)/(1-leaf.φ));

    leaf.Fo  = Kf/(Kf+Kp_max+Kd        );
    leaf.Fo′ = Kf/(Kf+Kp_max+Kd+leaf.Kn);
    leaf.Fm  = Kf/(Kf       +Kd        );
    leaf.Fm′ = Kf/(Kf       +Kd+leaf.Kn);
    leaf.ϕs  = leaf.Fm′*(1-leaf.φ);
    # leaf.eta  = leaf.ϕs/leaf.Fo
    # don't need this anymore
    # better to use ϕs directly for SIF as Fo is not always fqe=0.01
    leaf.qQ  = 1-(leaf.ϕs-leaf.Fo′)/(leaf.Fm-leaf.Fo′);
    leaf.qE  = 1-(leaf.Fm-leaf.Fo′)/(leaf.Fm′-leaf.Fo);
    leaf.NPQ = leaf.Kn/(Kf+Kd);

    return nothing
end
