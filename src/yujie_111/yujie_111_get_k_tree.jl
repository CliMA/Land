function Yujie111GetKTree(node::Yujie111{FT}) where {FT<:AbstractFloat}
    # 1. the root, stem, and leaf
    k_root = node.k_root
    k_stem = node.k_stem
    k_leaf = node.k_leaf
    l_root = node.l_root[:,2]
    l_stem = node.l_stem[:,2]
    l_leaf = node.l_leaf[:,2]
    r_root = 1.0/k_root * sum(1.0 ./ l_root) * 0.05
    r_stem = 1.0/k_stem * sum(1.0 ./ l_stem) * 0.05
    r_leaf = 1.0/k_leaf * sum(1.0 ./ l_leaf) * 0.05
    r_init = 1.0/k_root + 1.0/k_stem + 1.0/k_leaf

    # calculate the ks
    k_roo = 1.0/k_root / r_root
    k_ste = 1.0/k_stem / r_stem
    k_lea = 1.0/k_leaf / r_leaf
    k_all = r_init / (r_root + r_stem + r_leaf)
    return [k_all k_roo k_ste k_lea]
end
