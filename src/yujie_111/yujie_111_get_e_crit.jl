function Yujie111GetECrit(node::Yujie111{FT}) where {FT<:AbstractFloat}
    e_min = FT(0)
    e_max = FT(100)
    e_mid = FT(50)
    p_crt = -xylem_p_crit(node.vc_leaf);
    count = 0
    # make e_max 2X if not reaching p_crt
    while true
        count += 1
        if count>500
            break
        end
        p_max = Yujie111GetP(node, e_max)
        if p_max > p_crt
            break
        else
            e_min = e_max
            e_max *= 2
        end
    end
    # use bi-section method
    while true
        count += 1
        if count>500
            break
        end
        e_mid = 0.5 * (e_min + e_max)
        p_mid = Yujie111GetP(node, e_mid)
        if (abs(p_mid-p_crt)<1E-4) || (e_max-e_min<1E-4)
            break
        elseif p_mid > p_crt
            e_max = e_mid
        else
            e_min = e_mid
        end
    end
    return e_mid
end
