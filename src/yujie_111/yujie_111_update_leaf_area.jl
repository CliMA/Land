function Yujie111UpdateLeafArea(node::SPACSimple{FT}, laba) where {FT<:AbstractFloat}
    node.laba         = laba;
    node.hs.leaf.area = node.hs.stem.area * laba;
    node.lai          = laba / node.gaba;

    return nothing
end
