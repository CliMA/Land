
"""
    optimize_flows!(node::SPACSimple{FT}, photo_set::AbstractPhotoModelParaSet{FT}) where {FT<:AbstractFloat}

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
    node.ec = tree_e_crit(node.hs, node.ec);
    ec_sl   = node.ec * frac_sl;
    ec_sh   = node.ec * frac_sh;

    # optimize the sunlit and shade layers
    f_sl = min(ec_sl / FT(1.01), node.opt_f_sl);
    f_sh = min(ec_sh / FT(1.01), node.opt_f_sh);

    leaf_gas_exchange_nonopt!(node, photo_set, f_sl, f_sh)
    p_opt = node.containerOP

    # until df < 0.001
    df       = FT(0.1);
    f_tmp    = FT(0);
    p_tmp    = FT(0);
    count_sl = 0;
    count_sh = 0;
    while df>9e-4
        # increase the f_sl by df
        count_sl = 0
        while true
            f_tmp = f_sl + df;
            f_tmp > ec_sl ? break : nothing;
            leaf_gas_exchange_nonopt!(node, photo_set, f_tmp, f_sh);
            p_tmp = node.containerOP;
            p_tmp > p_opt ? (count_sl+=1; f_sl=f_tmp; p_opt=p_tmp;) : break;
        end

        # decrease the f_sl by df
        if count_sl == 0
            while true
                f_tmp = f_sl - df;
                f_tmp < 0 ? break : nothing;
                leaf_gas_exchange_nonopt!(node, photo_set, f_tmp, f_sh);
                p_tmp = node.containerOP;
                p_tmp > p_opt ? (count_sl+=1; f_sl=f_tmp; p_opt=p_tmp;) : break;
            end
        end

        # increase f_sh by df
        count_sh = 0
        while true
            f_tmp = f_sh + df;
            f_tmp > ec_sh ? break : nothing;
            leaf_gas_exchange_nonopt!(node, photo_set, f_sl, f_tmp);
            p_tmp = node.containerOP;
            p_tmp > p_opt ? (count_sh+=1; f_sh=f_tmp; p_opt=p_tmp;) : break;
        end

        # decrease f_sh by df
        if count_sh == 0
            while true
                f_tmp = f_sh - df;
                f_tmp < 0 ? break : nothing;
                leaf_gas_exchange_nonopt!(node, photo_set, f_sl, f_tmp);
                p_tmp = node.containerOP;
                p_tmp > p_opt ? (count_sh+=1; f_sh=f_tmp; p_opt=p_tmp;) : break;
            end
        end

        # 10% the df
        (count_sl==0 && count_sh==0) ? df /= 10 : nothing;
    end

    # return the flows
    node.opt_f_sl = f_sl;
    node.opt_f_sh = f_sh;

    return nothing
end
