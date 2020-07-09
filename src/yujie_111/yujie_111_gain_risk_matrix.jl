function Yujie111GainRiskMatrix(
            node::Yujie111{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            envir::AirLayer{FT},
            zenith=30.0,
            r_all=1000.0
) where {FT<:AbstractFloat}
    # partition the canopy
    canopy = Yujie111GetLeafPartition(node, zenith, r_all)
    r_sl = canopy[1]
    r_sh = canopy[4]

    # calculate the ecrit
    e_crit = Yujie111GetECrit(node)

    # mapping the matrix
    f_sl      = 0.0
    f_sh      = 0.0
    judge_sl  = true
    judge_sh  = true
    list_f_sl = 0.0
    list_f_sh = 0.0
    # initialize the matrix
    tmp    = Yujie111GetPACGTs(node, photo_set, f_sl, f_sh, canopy, envir)
    a_sum  = r_sl * tmp[1,3] + r_sh * tmp[2,3]
    e_sum  = f_sl + f_sh
    p_tmp  = (e_crit-e_sum) * a_sum
    matrix = p_tmp
    # expand the matrix
    while judge_sl || judge_sh
        # print f_sl and f_sh for debug
        #println("Sunlit: ", (@sprintf "%8.1f" f_sl), "\tShade: ", (@sprintf "%8.1f" f_sh))

        # if sunlit f can be higher
        tmp_list = []
        if judge_sl
            f_sl += 1.0
            # iterate through the shade layer f
            for tmp_f in list_f_sh
                tmp    = Yujie111GetPACGTs(node, photo_set, f_sl, tmp_f, canopy, envir)
                a_sum  = r_sl * tmp[1,3] + r_sh * tmp[2,3]
                e_sum  = f_sl + tmp_f
                p_tmp  = (e_crit-e_sum) * a_sum
                # judge if break
                if (tmp[1,3]<-1111.0)
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
            f_sh += 1.0
            # iterate through the sunlit layer f
            for tmp_f in list_f_sl
                tmp    = Yujie111GetPACGTs(node, photo_set, tmp_f, f_sh, canopy, envir)
                a_sum  = r_sl * tmp[1,3] + r_sh * tmp[2,3]
                e_sum  = tmp_f + f_sh
                p_tmp  = (e_crit-e_sum) * a_sum
                # judge if break
                if (tmp[2,3]<-1111.0)
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
