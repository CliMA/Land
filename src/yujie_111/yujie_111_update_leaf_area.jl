function Yujie111UpdateLeafArea(node::Yujie111{FT}, laba) where {FT<:AbstractFloat}
    node.laba   = laba
    node.k_leaf = node.k_sla * laba
    node.gmax   = laba * node.g_sla * 3600.0 * 18.0 * 0.001
end
