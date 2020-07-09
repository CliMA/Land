function Yujie111UpdateSoil(node::Yujie111{FT}, flow, dt=1.0) where {FT<:AbstractFloat}
    # 1. positive flow means out, flow in Kg h-1, dt in h
    m_ini = node.c_curr * node.gaba * node.h_soil * 998.0 # in Kg
    m_out = flow * dt
    m_end = m_ini - m_out
    c_end = m_end / (node.gaba * node.h_soil * 998.0)

    # 2. id c_end <0
    if c_end < 1E-6
        c_end = 1E-6
    end
    Yujie111UpdateSoilFromSWC(node, c_end)
end
