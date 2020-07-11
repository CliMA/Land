# this function is meant to evluate the models thru minor changes in flow
function Yujie111EvaluateModel(
            node::SPACSimple{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            zenith::FT,
            r_all::FT
) where {FT<:AbstractFloat}
    # partition the canopy
    big_leaf_partition!(node, zenith, r_all);

    # calculate the ecrit
    node.ec = tree_e_crit(node.hs, node.ec);
    p_crit = node.hs.leaf.p_crt

    # increase the flow
    df = DataFrame(E  = FT[],
                   P  = FT[],
                   An = FT[],
                   Ag = FT[],
                   C  = FT[],
                   G  = FT[],
                   T  = FT[],
                   Ec = FT[],
                   Pc = FT[])
    flow = FT(0);
    while flow<node.ec
        Yujie111GetPACGT(node, photo_set, flow)
        if (node.container1L).p < p_crit
            break
        end
        # expand the DataFrame
        push!(df, [ (node.container1L).e,
                    (node.container1L).p,
                    (node.container1L).an,
                    (node.container1L).ag,
                    (node.container1L).c,
                    (node.container1L).gh,
                    (node.container1L).t,
                    node.ec,
                    p_crit])
        flow += 0.1
        #println(flow, "\t", e_crit)
    end

    # return the matrix
    return df
end