"""
    get_relative_surface_tension(tem)
This function returns the surface tension of water relative to 25 degree C (298.15 K).
γ     = V^(-2/3) * k * (tem_c-tem)
γ/γ25 = (tem_c - tem) / (tem_c - 298.15)
The empirical values are
k     = 2.1E-7 J K^-1 mol^(-2/3)
V     = 18.0 ml/mol
tem_c = 647.0 K
"""
function get_relative_surface_tension(tem::FT) where {FT}
    return (st_tc - tem) / (st_tc - K_25)
end
