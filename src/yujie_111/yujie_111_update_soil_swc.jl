function Yujie111UpdateSoilFromSWC(node::SPACSimple{FT}, swc) where {FT<:AbstractFloat}
    if swc<node.c_ssat
        node.c_curr = swc
    else
        node.c_curr = node.c_ssat
    end

    node.p_soil = soil_p_25_swc(node.hs.root.sh, node.c_curr)
end
