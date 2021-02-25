###############################################################################
#
# Optimize leaf invstment Vcmax+Jmax and Leaf area
#
###############################################################################
"""
    optimize_hs!(
                node::SPACSimple{FT},
                photo_set::AbstractPhotoModelParaSet{FT},
                weather::Array{FT,2}
    ) where {FT<:AbstractFloat}

Optimize hydraulic conductance and leaf investment, given
- `node` [`SPACSimple`] type struct
- `photo_set` [`AbstractPhotoModelParaSet`] type struct
- `weather` Weather profile in a growing season
"""
function optimize_hs!(
            node::SPACSimple{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            weather::Array{FT,2}
) where {FT<:AbstractFloat}
    # 1. use the opt_laba and opt_vmax to initialize
    @inline f(x) = (tmp_node = deepcopy(node);
                    tmp_node.hs.root.k_max = tmp_node.containerKS[1] * x[1];
                    tmp_node.hs.stem.k_max = tmp_node.containerKS[2] * x[1];
                    tmp_node.hs.leaf.k_sla = tmp_node.containerKS[3] * x[1];
                    tmp_node.hs.root.k_element .= tmp_node.hs.root.k_max *
                                                  tmp_node.hs.root.N;
                    tmp_node.hs.stem.k_element .= tmp_node.hs.stem.k_max *
                                                  tmp_node.hs.stem.N;
                    tmp_node.hs.leaf.k_element .= tmp_node.hs.leaf.k_sla *
                                                  tmp_node.hs.leaf.N;
                    leaf_allocation!(tmp_node, photo_set, x[2], x[3]);
                    tmp_prof = annual_profit(tmp_node, photo_set, weather);
                    return tmp_prof);

    ms = ReduceStepMethodND{FT}(
                x_mins=FT[0.001,1,5],
                x_maxs=FT[1e10, node.gaba*10, 200],
                x_inis=FT[node.hs.root.k_max/node.containerKS[1],
                          node.opt_laba,
                          node.opt_vmax],
                Δ_inis=FT[0.1,100,10]);
    st = SolutionToleranceND{FT}(FT[0.0009, 0.9, 0.09], 50);
    klv = find_peak(f, ms, st);

    # 2. update the opt_laba and opt_vmax to node
    node.hs.root.k_max = node.containerKS[1] * klv[1];
    node.hs.stem.k_max = node.containerKS[2] * klv[1];
    node.hs.leaf.k_sla = node.containerKS[3] * klv[1];
    node.hs.root.k_element .= node.hs.root.k_max * node.hs.root.N;
    node.hs.stem.k_element .= node.hs.stem.k_max * node.hs.stem.N;
    node.hs.leaf.k_element .= node.hs.leaf.k_sla * node.hs.leaf.N;
    leaf_allocation!(node, photo_set, klv[2], klv[3]);
    node.opt_laba = klv[2];
    node.opt_vmax = klv[3];

    return nothing
end
