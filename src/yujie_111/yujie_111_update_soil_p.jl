function Yujie111UpdateSoilFromP(node::SPACSimple{FT}, p::FT) where {FT<:AbstractFloat}
    node.p_soil = p;
    node.c_curr = soil_swc(node.hs.root.sh, p);
end
