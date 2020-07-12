###############################################################################
#
# Calculate optimal flow rates into sunlit and shaded layers
#
###############################################################################

function flow_optimizer(node::SPACSimple{FT}, photo_set::AbstractPhotoModelParaSet{FT}, fs::Array{FT,1}) where {FT<:AbstractFloat}
    @unpack frac_sh, frac_sl = node.container2L;

    f_sl, f_sh = fs
    leaf_gas_exchange!(node, photo_set, f_sl, f_sh)
    a_sum = frac_sl * (node.container2L).cont_sl.an + frac_sh * (node.container2L).cont_sh.an;
    e_sum = f_sl + f_sh
    p_opt = (node.ec-e_sum) * a_sum

    return -p_opt
end




function optimize_flows_testing(
            node::SPACSimple{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            zenith::FT,
            r_all::FT
) where {FT<:AbstractFloat}
    # partition the canopy
    big_leaf_partition!(node, zenith, r_all)
    @unpack frac_sh, frac_sl = node.container2L;

    # calculate the ecrit
    node.ec = tree_e_crit(node.hs, node.ec);

    # optimize the sunlit and shade layers using initial values
    f_sl = frac_sl * node.ec / 3
    f_sh = frac_sh * node.ec / 5

    @inline f(x) = flow_optimizer(node, photo_set, x)
    _res = optimize(f, [f_sl, f_sh], NelderMead())

    f_sl = Optim.minimizer(_res)[1];
    f_sh = Optim.minimizer(_res)[2];

    # return the flows
    node.opt_f_sl = f_sl;
    node.opt_f_sh = f_sh;

    return nothing
end