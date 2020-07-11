function Yujie111UpdateLeaf(node::SPACSimple{FT}, photo_set::AbstractPhotoModelParaSet{FT}, laba, vmax, jmaxr=1.75) where {FT<:AbstractFloat}
    node.laba         = laba;
    node.hs.leaf.area = node.hs.stem.area * laba;
    node.lai          = laba / node.gaba;
    node.ps.Vcmax25   = vmax;
    node.ps.Jmax25    = vmax * jmaxr;
    node.ps.Rd25      = vmax * photo_set.VR;

    return nothing
end
