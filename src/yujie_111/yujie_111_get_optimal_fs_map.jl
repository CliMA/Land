# TODO remove the all -Inf lines




function Yujie111GetOptimalFsMap(
            node::SPACSimple{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            zenith=30.0,
            r_all=1000.0
) where {FT<:AbstractFloat}
    # partition the canopy
    big_leaf_partition!(node, zenith, r_all)
    @unpack frac_sh, frac_sl = node.container2L;

    # calculate the ecrit
    node.ec = tree_e_crit(node.hs, node.ec);

    # make matrix for plotting
    mat  = []
    f_sl = 0.0
    while true
        #println(f_sl, " in ", node.ec)
        f_sh = 0.0
        row  = []
        while true
            Yujie111GetPACGTs(node, photo_set, f_sl, f_sh)
            a_sum = frac_sl * node.container2L.cont_sl.an + frac_sh * node.container2L.cont_sh.an;
            e_sum = f_sl + f_sh
            p_opt = (node.ec-e_sum) * a_sum / node.ec
            if row==[]
                row = p_opt
            else
                row = [row; p_opt]
            end
            f_sh += 0.1
            if f_sh>node.ec
                break
            end
        end
        if mat==[]
            mat = row
        else
            mat = [mat row]
        end
        f_sl += 0.1
        if f_sl>node.ec
            break
        end
    end

    # return the map
    return mat
end
