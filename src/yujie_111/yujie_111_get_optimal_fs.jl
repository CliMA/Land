function Yujie111GetOptimalFs(
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

    # optimize the sunlit and shade layers
    f_sl  = 0.0
    f_sh  = 0.0
    tmp   = Yujie111GetPACGTs(node, photo_set, f_sl, f_sh, canopy, envir)
    a_sum = r_sl * tmp[1,3] + r_sh * tmp[2,3]
    e_sum = f_sl + f_sh
    p_opt = (e_crit-e_sum) * a_sum
    df = 100.0
    while df>0.1
        # print the f_sl and f_sh and df
        #println( (@sprintf "df: %8.3f\tf_sl: %8.3f\tf_sh: %8.3f" df f_sl f_sh) )

        # increase the f_sl by df
        count = 0
        if count<1
            while true
                f_tmp = max(0.0, f_sl+df)
                if f_tmp+f_sh>e_crit
                    break
                end
                tmp   = Yujie111GetPACGTs(node, photo_set, f_tmp, f_sh, canopy, envir)
                a_sum = r_sl * tmp[1,3] + r_sh * tmp[2,3]
                e_sum = f_tmp + f_sh
                p_tmp = (e_crit-e_sum) * a_sum
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
                if f_tmp+f_sl>e_crit
                    break
                end
                tmp   = Yujie111GetPACGTs(node, photo_set, f_tmp, f_sh, canopy, envir)
                a_sum = r_sl * tmp[1,3] + r_sh * tmp[2,3]
                e_sum = f_tmp + f_sh
                p_tmp = (e_crit-e_sum) * a_sum
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
                if f_sl+f_tmp>e_crit
                    break
                end
                tmp   = Yujie111GetPACGTs(node, photo_set, f_sl, f_tmp, canopy, envir)
                a_sum = r_sl * tmp[1,3] + r_sh * tmp[2,3]
                e_sum = f_sl + f_tmp
                p_tmp = (e_crit-e_sum) * a_sum
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
                if f_sl+f_tmp>e_crit
                    break
                end
                tmp   = Yujie111GetPACGTs(node, photo_set, f_sl, f_tmp, canopy, envir)
                a_sum = r_sl * tmp[1,3] + r_sh * tmp[2,3]
                e_sum = f_sl + f_tmp
                p_tmp = (e_crit-e_sum) * a_sum
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
    return [f_sl f_sh]
end
