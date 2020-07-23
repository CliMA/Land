###############################################################################
#
# Optimize leaf invstment Vcmax+Jmax and Leaf area
#
###############################################################################
"""
    optimize_leaf!(node::SPACSimple{FT}, photo_set::AbstractPhotoModelParaSet{FT}, weather::Array{FT,2}) where {FT<:AbstractFloat}

Optimize leaf area and photosynthetic capacity, given
- `node` [`SPACSimple`] type struct
- `photo_set` [`AbstractPhotoModelParaSet`] type struct
- `weather` Weather profile in a growing season
"""
function optimize_leaf!(
            node::SPACSimple{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            weather::Array{FT,2}
) where {FT<:AbstractFloat}
    # 1. use the opt_laba and opt_vmax to initialize
    @inline f(x) = (tmp_node = deepcopy(node);
                    leaf_allocation!(tmp_node, photo_set, x[1], x[2]);
                    tmp_prof = annual_profit(tmp_node, photo_set, weather);
                    return tmp_prof);

    ms = ReduceStepMethodND{FT}(FT[0,0],
                                FT[1e10, node.maxv],
                                FT[node.opt_laba, node.opt_vmax],
                                FT[100,10]);
    st = SolutionToleranceND{FT}(FT[0.9, 0.09], 50);
    lv = find_peak(f, ms, st);

    # 2. update the opt_laba and opt_vmax to node
    leaf_allocation!(node, photo_set, lv[1], lv[2]);
    node.opt_laba = lv[1];
    node.opt_vmax = lv[2];

    return nothing
end
