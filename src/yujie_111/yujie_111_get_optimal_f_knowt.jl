# function to test birch dataset, know t_leaf
function Yujie111GetOptimalFKnowT(
            node::Yujie111{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            envir::AirLayer{FT},
            conditions
) where {FT<:AbstractFloat}
    # get e_crit
    e_crit = Yujie111GetECrit(node)

    # use bi-section method
    e_min = 0.0
    e_max = e_crit
    r_opt = []
    while true
        e_mid       = 0.5 * (e_min+e_max)
        p,a,c,g     = Yujie111GetPACGKnowT(node, photo_set, envir, conditions, e_mid)
        e_de        = e_mid + 1.0
        pd,ad,cd,gd = Yujie111GetPACGKnowT(node, photo_set, envir, conditions, e_de )
        optimizer   = ad*(e_crit-e_de) - a*(e_crit-e_mid)
        if (e_max-e_min)<1.0




            r_opt = [e_mid p a c g conditions[4]]




            break
        end
        if optimizer>0
            e_min = e_mid
        else
            e_max = e_mid
        end
    end

    # return the optima
    return r_opt
end
