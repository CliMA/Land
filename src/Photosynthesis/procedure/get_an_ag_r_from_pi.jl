"""
    get_an_ag_r_from_pi(model::C3ParaSetVcJ, p_i::FT, v25::FT, j25::FT, par::FT, p_O₂::FT, t_leaf::FT, r25::FT, curvature::FT, qy::FT)

Gross photosynthetic rate limited by carboxylation without ATP limitation, given
- `model` A C3ParaSetVcJ type parameter sets that stores temperature denpendencies
- `p_i` Leaf internal CO₂ partial pressure
- `v25` Maximal carboxylation rate at 298.15 K
- `j25` Maximal electron transport rate
- `par` Absorbed photosynthetic active radiation
- `p_O₂` O₂ partial pressure
- `t_leaf` Leaf temperature
- `r25` Respiration rate at 298.15 K
- `curvature` Curvature parameter to calculate J
- `qy` Quantum yield of electron

The equations used are from Farquhar et al. (1980) "A biochemical model of photosynthetic CO2 assimilation in leaves of C3 species."

"""
function get_an_ag_r_from_pi(model::C3ParaSetVcJ,
                               p_i::FT,
                               v25::FT,
                               j25::FT,
                               par::FT,
                              p_O₂::FT,
                            t_leaf::FT,
                               r25::FT,
                         curvature::FT,
                                qy::FT) where {FT}
    if r25==Inf
        r25 = model.VR.v_to_r * v25
    end
    r_leaf  = get_r(model.ReT, r25, t_leaf)
    vmax    = get_vmax(model.VcT, v25, t_leaf)
    jmax    = get_jmax(model.JT , j25, t_leaf)
    j       = get_j(jmax, par, curvature, qy)
    kc      = get_kc(model.KcT, t_leaf)
    ko      = get_ko(model.KoT, t_leaf)
    km      = kc * (1 + p_O₂/ko)
    Γ_star  = get_Γ_star(model.ΓsT, t_leaf)
    a_v     = vmax * (p_i-Γ_star) / (p_i+km)
    a_j     = j * (p_i-Γ_star) / ( 4*(p_i + 2*Γ_star) )
    ag      = min(a_v, a_j)
    an      = ag - r_leaf
    return an,ag,r_leaf
end




"""
    get_an_ag_r_from_pi(model::C3ParaSetVcVpJ, p_i::FT, v25::FT, j25::FT, par::FT, p_O₂::FT, t_leaf::FT, r25::FT, curvature::FT, qy::FT)

Gross photosynthetic rate limited by carboxylation with ATP limitation, given
- `model` A C3ParaSetVcVpJ type parameter sets that stores temperature denpendencies
- `p_i` Leaf internal CO₂ partial pressure
- `v25` Maximal carboxylation rate at 298.15 K
- `j25` Maximal electron transport rate
- `par` Absorbed photosynthetic active radiation
- `p_O₂` O₂ partial pressure
- `t_leaf` Leaf temperature
- `r25` Respiration rate at 298.15 K
- `curvature` Curvature parameter to calculate J
- `qy` Quantum yield of electron

The equations used are from Farquhar et al. (1980) "A biochemical model of photosynthetic CO2 assimilation in leaves of C3 species."

"""
function get_an_ag_r_from_pi(model::C3ParaSetVcVpJ,
                               p_i::FT,
                               v25::FT,
                               j25::FT,
                               par::FT,
                              p_O₂::FT,
                            t_leaf::FT,
                               r25::FT,
                         curvature::FT,
                                qy::FT) where {FT}
    if r25==Inf
        r25 = model.VR.v_to_r * v25
    end
    r_leaf  = get_r(model.ReT, r25, t_leaf)
    vmax    = get_vmax(model.VcT, v25, t_leaf)
    jmax    = get_jmax(model.JT, j25, t_leaf)
    j       = get_j(jmax, par, curvature, qy)
    kc      = get_kc(model.KcT, t_leaf)
    ko      = get_ko(model.KoT, t_leaf)
    km      = kc * (1 + p_O₂/ko)
    Γ_star  = get_Γ_star(model.ΓsT, t_leaf)
    a_p     = vmax / 2
    a_v     = vmax * (p_i-Γ_star) / (p_i+km)
    a_j     = j * (p_i-Γ_star) / ( 4*(p_i + 2*Γ_star) )
    ag      = min(a_p, a_v, a_j)
    an      = ag - r_leaf
    return an,ag,r_leaf
end




"""
    get_an_ag_r_from_pi(model::C4ParaSetVcVpJ, p_i::FT, v25::FT, p25::FT, par::FT, t_leaf::FT, r25::FT, qy::FT)

Gross photosynthetic rate limited by carboxylation, given
- `model` A C4ParaSetVcVpJ type parameter sets that stores temperature denpendencies
- `p_i` Leaf internal CO₂ partial pressure
- `v25` Maximal carboxylation rate at 298.15 K
- `p25` Maximal PEP carboxylation rate at 298.15 K
- `par` Absorbed photosynthetic active radiation
- `t_leaf` Leaf temperature
- `r25` Respiration rate at 298.15 K
- `qy` Quantum yield of electron

The model is adapted from Collatz et al. (1992) "Coupled photosynthesis-stomatal conductance model for leaves of C4 plants."

"""
function get_an_ag_r_from_pi(model::C4ParaSetVcVpJ,
                               p_i::FT,
                               v25::FT,
                               p25::FT,
                               par::FT,
                            t_leaf::FT,
                               r25::FT,
                                qy::FT) where {FT}
    if r25==Inf
        r25 = model.VR.v_to_r * v25
    end
    r_leaf  = get_r(model.ReT, r25, t_leaf)
    pmax    = get_vmax(model.VpT, p25, t_leaf)
    a_v     = get_vmax(model.VcT, v25, t_leaf)
    kpep    = get_kpep(model.KpT, t_leaf)
    # different a_j calculation from the C3 model
    a_j     = qy * par / 6
    a_p     = pmax * p_i / (p_i + kpep)
    ag      = min(a_v, a_j, a_p)
    an      = ag - r_leaf
    return an,ag,r_leaf
end
