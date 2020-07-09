function Yujie111UpdateSoilFromP(node::Yujie111{FT}, p) where {FT<:AbstractFloat}
    if p>node.p_ssat
        node.p_soil = p
        node.c_curr = node.c_ssat * (p/node.p_ssat)^(-1.0/node.b_ssat)
    else
        node.p_soil = node.p_ssat
        node.c_curr = node.c_ssat
    end
end
