# Use more advanced algorithm?




function Yujie111GetOptimalFs(
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

    # optimize the sunlit and shade layers
    f_sl  = 0.0
    f_sh  = 0.0
    Yujie111GetPACGTs(node, photo_set, f_sl, f_sh)
    a_sum = frac_sl * node.container2L.cont_sl.an + frac_sh * node.container2L.cont_sh.an;
    e_sum = f_sl + f_sh
    p_opt = (node.ec-e_sum) * a_sum
    df = 2.0
    while df>0.0019
        # print the f_sl and f_sh and df
        #println( (@sprintf "df: %8.3f\tf_sl: %8.3f\tf_sh: %8.3f" df f_sl f_sh) )

        # increase the f_sl by df
        count = 0
        if count<1
            while true
                f_tmp = max(0.0, f_sl+df)
                if f_tmp+f_sh>node.ec
                    break
                end
                Yujie111GetPACGTs(node, photo_set, f_tmp, f_sh)
                a_sum = frac_sl * node.container2L.cont_sl.an + frac_sh * node.container2L.cont_sh.an;
                e_sum = f_tmp + f_sh
                p_tmp = (node.ec-e_sum) * a_sum
                if p_tmp>p_opt
                    #println("sl+", tmp)
                    count += 1
                    f_sl   = f_tmp
                    p_opt  = p_tmp
                else
                    break
                end
            end
        end

        # decrease the f_sl by df
        if count<1
            while true
                f_tmp = max(0.0, f_sl-df)
                if f_tmp+f_sl>node.ec
                    break
                end
                Yujie111GetPACGTs(node, photo_set, f_tmp, f_sh)
                a_sum = frac_sl * node.container2L.cont_sl.an + frac_sh * node.container2L.cont_sh.an;
                e_sum = f_tmp + f_sh
                p_tmp = (node.ec-e_sum) * a_sum
                if p_tmp>p_opt
                    #println("sl-", tmp)
                    count += 1
                    f_sl   = f_tmp
                    p_opt  = p_tmp
                else
                    break
                end
            end
        end

        # increase the f_sh by df
        count = 0
        if count<1
            while true
                f_tmp = max(0.0, f_sh+df)
                if f_sl+f_tmp>node.ec
                    break
                end
                Yujie111GetPACGTs(node, photo_set, f_sl, f_tmp)
                a_sum = frac_sl * node.container2L.cont_sl.an + frac_sh * node.container2L.cont_sh.an;
                e_sum = f_sl + f_tmp
                p_tmp = (node.ec-e_sum) * a_sum
                if p_tmp>p_opt
                    #println("sh+", tmp)
                    count += 1
                    f_sh   = f_tmp
                    p_opt  = p_tmp
                else
                    break
                end
            end
        end

        # decrease the f_sh by df
        if count<1
            while true
                f_tmp = max(0.0, f_sh-df)
                if f_sl+f_tmp>node.ec
                    break
                end
                Yujie111GetPACGTs(node, photo_set, f_sl, f_tmp)
                a_sum = frac_sl * node.container2L.cont_sl.an + frac_sh * node.container2L.cont_sh.an;
                e_sum = f_sl + f_tmp
                p_tmp = (node.ec-e_sum) * a_sum
                if p_tmp>p_opt
                    #println("sh-", tmp)
                    count += 1
                    f_sh   = f_tmp
                    p_opt  = p_tmp
                else
                    break
                end
            end
        end

        # half the df
        df *= 0.1
    end

    # return the flows
    node.opt_f_sl = f_sl;
    node.opt_f_sh = f_sh;

    return nothing
end
