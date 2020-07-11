function Yujie111GainRiskMatrix(
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

    # mapping the matrix
    f_sl      = 0.0
    f_sh      = 0.0
    judge_sl  = true
    judge_sh  = true
    list_f_sl = 0.0
    list_f_sh = 0.0

    # initialize the matrix
    Yujie111GetPACGTs(node, photo_set, f_sl, f_sh)
    a_sum  = frac_sl * node.container2L.cont_sl.an + frac_sh * node.container2L.cont_sh.an;
    e_sum  = f_sl + f_sh
    p_tmp  = (node.ec-e_sum) * a_sum / node.ec
    matrix = p_tmp

    # expand the matrix
    while judge_sl || judge_sh
        # print f_sl and f_sh for debug
        #println("Sunlit: ", (@sprintf "%8.1f" f_sl), "\tShade: ", (@sprintf "%8.1f" f_sh))

        # if sunlit f can be higher
        tmp_list = []
        if judge_sl
            f_sl += 0.1
            # iterate through the shade layer f
            for tmp_f in list_f_sh
                Yujie111GetPACGTs(node, photo_set, f_sl, tmp_f)
                a_sum  = frac_sl * node.container2L.cont_sl.an + frac_sh * node.container2L.cont_sh.an;
                e_sum  = f_sl + tmp_f
                p_tmp  = (node.ec-e_sum) * a_sum / node.ec

                # judge if break
                if (node.container2L.cont_sl.an < -1111.0)
                    judge_sl = false
                    break
                end
                # expand the tmp_list
                if tmp_list==[]
                    tmp_list = p_tmp
                else
                    tmp_list = [tmp_list p_tmp]
                end
            end
        end
        # if expand the matrix
        if judge_sl
            matrix    = [matrix; tmp_list]
            list_f_sl = [list_f_sl f_sl]
        end

        # if shade f can be higher
        tmp_list = []
        if judge_sh
            f_sh += 0.1
            # iterate through the sunlit layer f
            for tmp_f in list_f_sl
                Yujie111GetPACGTs(node, photo_set, tmp_f, f_sh)
                a_sum  = frac_sl * node.container2L.cont_sl.an + frac_sh * node.container2L.cont_sh.an;
                e_sum  = tmp_f + f_sh
                p_tmp  = (node.ec-e_sum) * a_sum / node.ec

                # judge if break
                if (node.container2L.cont_sh.an < -1111.0)
                    judge_sh = false
                    break
                end
                # expand the tmp_list
                if tmp_list==[]
                    tmp_list = p_tmp
                else
                    tmp_list = [tmp_list; p_tmp]
                end
            end
        end
        # if expand the matrix
        if judge_sh
            matrix    = [matrix tmp_list]
            list_f_sh = [list_f_sh f_sh]
        end
    end

    # return the matrix
    return matrix
end
