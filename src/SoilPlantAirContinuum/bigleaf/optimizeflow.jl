
###############################################################################
#
# Optimize flow rates in the sunlit and shaded layers
#
###############################################################################
"""
    optimize_flows!(
                node::SPACSimple{FT},
                photo_set::AbstractPhotoModelParaSet{FT}
    ) where {FT<:AbstractFloat}

Optimize the flow rates in sunlit and shaded layers, given
- `node` [`SPACSimple`] type struct
- `photo_set` [`AbstractPhotoModelParaSet`] type struct
"""
function optimize_flows!(
            node::SPACSimple{FT},
            photo_set::AbstractPhotoModelParaSet{FT}
) where {FT<:AbstractFloat}
    # unpack required parameters
    @unpack frac_sh, frac_sl = node.container2L;

    # calculate the ecrit
    node.ec = critical_flow(node.hs, node.ec);
    ec_sl   = node.ec * frac_sl;
    ec_sh   = node.ec * frac_sh;

    # optimize the sunlit and shade layers
    f_sl = min(ec_sl / FT(1.01), node.opt_f_sl);
    f_sh = min(ec_sh / FT(1.01), node.opt_f_sh);

    ms = ReduceStepMethodND{FT}(
                x_mins=FT[0,0],
                x_maxs=[ec_sl,ec_sh],
                x_inis=[f_sl, f_sh],
                Δ_inis=FT[0.1,0.1]);
    st = SolutionToleranceND{FT}(FT[9e-4, 9e-4], 50);
    @inline f(x) = (leaf_gas_exchange_nonopt!(node, photo_set, x[1], x[2]);
                    return node.containerOP);
    fs = find_peak(f, ms, st);

    # update the optimal flow rates
    node.opt_f_sl = fs[1];
    node.opt_f_sh = fs[2];

    return nothing
end
