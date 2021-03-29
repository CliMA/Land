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
                photo_set::AbstractPhotoModelParaSet{FT},
                vmax::FT
    ) where {FT<:AbstractFloat}
    leaf_allocation!(
                node::SPACSimple{FT},
                photo_set::AbstractPhotoModelParaSet{FT},
                laba::FT,
                vmax::FT
    ) where {FT<:AbstractFloat}

Update leaf area and maximal carboxylation rate, given
- `node` [`SPACSimple`] type struct
- `photo_set` [`AbstractPhotoModelParaSet`] type struct
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
            photo_set::AbstractPhotoModelParaSet{FT},
            vmax::FT
) where {FT<:AbstractFloat}
    node.ps.Vcmax25 = vmax;
    node.ps.Jmax25  = vmax * node.vtoj;
    node.ps.Rd25    = vmax * photo_set.VR;

    return nothing
end




function leaf_allocation!(
            node::SPACSimple{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            laba::FT,
            vmax::FT
) where {FT<:AbstractFloat}
    leaf_allocation!(node, laba);
    leaf_allocation!(node, photo_set, vmax);

    return nothing
end
