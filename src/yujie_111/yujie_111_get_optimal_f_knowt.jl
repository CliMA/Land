# function to test birch dataset, know t_leaf
# not yet tested




function Yujie111GetOptimalFKnowT(
            node::SPACSimple{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            conditions
) where {FT<:AbstractFloat}
    # get e_crit
    node.ec = tree_e_crit(node, node.ec);

    # use bi-section method
    e_min = 0.0
    e_max = node.ec
    r_opt = []
    while true
        e_mid       = 0.5 * (e_min+e_max)
        p,a,c,g     = Yujie111GetPACGKnowT(node, photo_set, conditions, e_mid)
        e_de        = e_mid + 1.0
        pd,ad,cd,gd = Yujie111GetPACGKnowT(node, photo_set, conditions, e_de )
        optimizer   = ad*(node.ec-e_de) - a*(node.ec-e_mid)
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
