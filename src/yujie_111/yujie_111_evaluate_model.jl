# this function is meant to evluate the models thru minor changes in flow
function Yujie111EvaluateModel(
            node::Yujie111{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            envir::AirLayer{FT},
            zenith::FT=FT(30),
            r_all::FT=FT(1000)
) where {FT<:AbstractFloat}
    # partition the canopy
    lai    = node.laba / node.gaba;
    canopy = zeros(FT,6);
    big_leaf_partition!(canopy, zenith, lai, r_all);

    r_sl = canopy[1]
    r_sh = canopy[4]
    @show canopy;

    # calculate the ecrit
    e_crit = Yujie111GetECrit(node);
    @show e_crit;
    p_crit = -xylem_p_crit(node.vc_leaf);
    @show p_crit;

    # increase the flow
    matrix = []
    flow = FT(0);
    while flow<e_crit
        @show flow;
        e,p,a,c,g,t = Yujie111GetPACGT(node, photo_set, flow, canopy, envir)
        if p<0
            break
        end
        # expand the matrix
        if matrix==[]
            matrix = [e p a c g t e_crit p_crit]
        else
            matrix = [matrix; e p a c g t e_crit p_crit]
        end
        flow += 0.1
        println( (@sprintf "%8.3f  in  %8.3f" flow e_crit) )
    end

    # return the matrix
    return matrix
end