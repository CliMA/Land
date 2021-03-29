###############################################################################
#
# Create a matrix to show the sunlit-shaded partition algorithm
#
###############################################################################
"""
    gain_risk_map(
                node::SPACSimple{FT},
                photo_set::AbstractPhotoModelParaSet{FT}
    ) where {FT<:AbstractFloat}

Return the matrix of optimizers at different sunlit and shaded layer flow
    rates, given
- `node` [`SPACSimple`] type struct
- `photo_set` [`AbstractPhotoModelParaSet`] type photosynthesis parameter set
"""
function gain_risk_map(
            node::SPACSimple{FT},
            photo_set::AbstractPhotoModelParaSet{FT}
) where {FT<:AbstractFloat}
    # 0. unpack required variables
    @unpack frac_sh, frac_sl = node.container2L;

    # 1. calculate the critical flow rate
    node.ec = critical_flow(node.hs, node.ec);

    # 2. find the maximal rates
    max_f_sl = node.ec * frac_sl;
    max_f_sh = node.ec * frac_sh;
    f_sl = FT(0);
    f_sh = FT(0);
    list_f_sl = FT[0];
    list_f_sh = FT[0];

    # 3. initialize the matrix
    leaf_gas_exchange_nonopt!(node, photo_set, f_sl, f_sh);
    mat = node.containerOP /  node.ec;

    # 4. expand the matrix
    judge_sl = true;
    judge_sh = true;
    Δf = FT(0.1);
    while judge_sl || judge_sh
        # 4.1 if sunlit f_sl can be higher
        tmp_list = [];
        if judge_sl
            f_sl += Δf;
            f_sl >= max_f_sl ? (judge_sl=false; break) : nothing;

            # 4.1.1 iterate through the list_f_sh
            for tmp_f_sh in list_f_sh
                leaf_gas_exchange_nonopt!(node, photo_set, f_sl, tmp_f_sh);
                ele = node.containerOP / node.ec;

                # expand the tmp_list horizontally
                if tmp_list==[]
                    tmp_list = ele;
                else
                    tmp_list = [tmp_list ele];
                end
            end

            # 4.1.2 break if all -Inf
            if all( tmp_list .< -999 )
                judge_sl = false;
            end
        end

        # 4.2 expand the matrix vertically
        if judge_sl
            mat = [mat; tmp_list];
            push!(list_f_sl, f_sl);
        end

        # 4.3 if shaded f_sh can be higher
        tmp_list = []
        if judge_sh
            f_sh += Δf;
            f_sh >= max_f_sh ? (judge_sh=false; break) : nothing;

            # 4.3.1 iterate through the list_f_sl
            for tmp_f_sl in list_f_sl
                leaf_gas_exchange_nonopt!(node, photo_set, tmp_f_sl, f_sh);
                ele = node.containerOP / node.ec;

                # expand the tmp_list vertically
                push!(tmp_list, ele);
            end

            # 4.3.2 break if all -Inf
            if all( tmp_list .< -999 )
                judge_sh = false;
            end
        end

        # 4.4 expand the matrix horizontally
        if judge_sh
            mat = [mat tmp_list];
            push!(list_f_sh, f_sh);
        end
    end

    # return the matrix, size(mat)[1] for sunlit, [2] for shaded
    return mat
end
