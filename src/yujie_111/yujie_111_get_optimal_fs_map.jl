function Yujie111GetOptimalFsMap(
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

    # make matrix for plotting
    mat  = []
    f_sl = 0.0
    while true
        println(f_sl, " in ", e_crit)
        f_sh = 0.0
        row  = []
        while true
            tmp   = Yujie111GetPACGTs(node, photo_set, f_sl, f_sh, canopy, envir)
            a_sum = r_sl * tmp[1,3] + r_sh * tmp[2,3]
            e_sum = f_sl + f_sh
            p_opt = (e_crit-e_sum) * a_sum / e_crit
            if row==[]
                row = p_opt
            else
                row = [row; p_opt]
            end
            f_sh += 1.0
            if f_sh>e_crit
                break
            end
        end
        if mat==[]
            mat = row
        else
            mat = [mat row]
        end
        f_sl += 1.0
        if f_sl>e_crit
            break
        end
    end

    # return the map
    return mat
end
