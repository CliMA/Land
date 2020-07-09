function Yujie111UpdateSoilFromSWC(node::Yujie111{FT}, swc) where {FT<:AbstractFloat}
    if swc<node.c_ssat
        node.c_curr = swc
        node.p_soil = node.p_ssat * (node.c_curr / node.c_ssat)^(-node.b_ssat)
    else
        node.c_curr = node.c_ssat
        node.p_soil = node.p_ssat
    end
end
