function Yujie111UpdateLeaf(node::Yujie111{FT}, laba, vmax, jmaxr=1.75) where {FT<:AbstractFloat}
    node.laba       = laba
    node.k_leaf     = node.k_sla * laba
    node.ps.Vcmax25 = vmax
    node.ps.Jmax25  = vmax * jmaxr
    node.gmax       = laba * node.g_sla * 3600.0 * 18.0 * 0.001
end
