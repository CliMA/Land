#=
###############################################################################
#
# Update leaf allocation
#
###############################################################################
"""
    leaf_allocation!(
                node::SPACSimple{FT},
                laba::FT
    ) where {FT<:AbstractFloat}
    leaf_allocation!(
                node::SPACSimple{FT},
                isvmax::Bool,
                vmax::FT,
                v2r::FT = FT(0.015)
    ) where {FT<:AbstractFloat}
    leaf_allocation!(
                node::SPACSimple{FT},
                laba::FT,
                vmax::FT
    ) where {FT<:AbstractFloat}

Update leaf area and maximal carboxylation rate, given
- `node` [`SPACSimple`] type struct
- `laba` Given leaf area per basal area
- `vmax` Given Vcmax25
"""
function leaf_allocation!(
            node::SPACSimple{FT},
            laba::FT
) where {FT<:AbstractFloat}
    node.laba         = laba;
    node.hs.leaf.area = node.hs.stem.area * laba;
    node.lai          = laba / node.gaba;

    return nothing
end




function leaf_allocation!(
            node::SPACSimple{FT},
            isvmax::Bool,
            vmax::FT,
            v2r::FT = FT(0.015)
) where {FT<:AbstractFloat}
    node.ps.PSM.v_cmax25_ww = vmax;
    node.ps.PSM.v_cmax25    = vmax;
    node.ps.PSM.j_max25     = vmax * node.vtoj;
    node.ps.PSM.r_d25       = vmax * v2r;

    return nothing
end




function leaf_allocation!(
            node::SPACSimple{FT},
            laba::FT,
            vmax::FT
) where {FT<:AbstractFloat}
    leaf_allocation!(node, laba);
    leaf_allocation!(node, true, vmax);

    return nothing
end
=#
