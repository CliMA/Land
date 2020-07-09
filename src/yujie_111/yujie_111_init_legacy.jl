function Yujie111InitializeLegacy(node::Yujie111{FT}) where {FT<:AbstractFloat}
    node.l_root  = [zeros(20) ones(20)]
    node.l_stem  = [zeros(20) ones(20)]
    node.l_leaf  = [zeros(20) ones(20)]
end
