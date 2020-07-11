function Yujie111UpdateSoil(node::SPACSimple{FT}, flow, dt=1.0) where {FT<:AbstractFloat}
    # 1. positive flow means out, flow in mol s-1, dt in h
    m_ini = node.c_curr * node.gaba * node.h_soil * ρ_H₂O # in Kg
    m_out = flow * dt
    m_end = m_ini - m_out
    c_end = m_end / (node.gaba * node.h_soil * ρ_H₂O)

    # 2. if c_end <0
    c_end = max(c_end, FT(1e-6))

    Yujie111UpdateSoilFromSWC(node, c_end)
    
    return nothing
end
