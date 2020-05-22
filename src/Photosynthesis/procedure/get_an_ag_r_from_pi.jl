"""
    get_an_ag_r_from_pi(; model, v25, j25, p_i, par, p_O₂, t_leaf, r25, curvature, qy)

Gross photosynthetic rate limited by carboxylation, given
- `model` A C3ParaSetVcJ type parameter sets that stores temperature denpendencies
- `v25` Maximal carboxylation rate at 298.15 K
- `j25` Maximal electron transport rate
- `p_i` Leaf internal CO₂ partial pressure
- `par` Absorbed photosynthetic active radiation
- `p_O₂` O₂ partial pressure
- `t_leaf` Leaf temperature
- `r25` Respiration rate at 298.15 K
- `curvature` Curvature parameter to calculate J
- `qy` Quantum yield of electron

The equations used are from Farquhar et al. (1980) "A biochemical model of photosynthetic CO2 assimilation in leaves of C3 species."

"""
function get_an_ag_r_from_pi(;
                             model::C3ParaSetVcJ = C3VcJBernacchi{FT}(),
                               v25::FT           = FT(80.0),
                               j25::FT           = FT(135.0),
                               p_i::FT           = FT(30.0),
                               par::FT           = FT(1000.0),
                              p_O₂::FT           = FT(21278.25),
                            t_leaf::FT           = FT(298.15),
                               r25::FT           = FT(Inf),
                         curvature::FT           = FT(0.9),
                                qy::FT           = FT(0.4081632653061224)) where {FT}
    if r25==Inf
        r25 = model.VR.v_to_r * v25
    end
    r_leaf  = get_r(model.ReT, r25, t_leaf)
    vmax    = get_vmax(model.VcT, v25, t_leaf)
    jmax    = get_jmax(model.JT, j25, t_leaf)
    j       = get_j(jmax, par; curvature=curvature, qy=qy)
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
    get_an_ag_r_from_pi(; model, v25, j25, p_i, par, p_O₂, t_leaf, r25, curvature, qy)

Gross photosynthetic rate limited by carboxylation, given
- `model` A C3ParaSetVcVpJ type parameter sets that stores temperature denpendencies
- `v25` Maximal carboxylation rate at 298.15 K
- `j25` Maximal electron transport rate
- `p_i` Leaf internal CO₂ partial pressure
- `par` Absorbed photosynthetic active radiation
- `p_O₂` O₂ partial pressure
- `t_leaf` Leaf temperature
- `r25` Respiration rate at 298.15 K
- `curvature` Curvature parameter to calculate J
- `qy` Quantum yield of electron

The equations used are from Farquhar et al. (1980) "A biochemical model of photosynthetic CO2 assimilation in leaves of C3 species."

"""
function get_an_ag_r_from_pi(;
                             model::C3ParaSetVcVpJ = C3VcJBernacchi{FT}(),
                               v25::FT             = FT(80.0),
                               j25::FT             = FT(135.0),
                               p_i::FT             = FT(30.0),
                               par::FT             = FT(1000.0),
                              p_O₂::FT             = FT(21278.25),
                            t_leaf::FT             = FT(298.15),
                               r25::FT             = FT(Inf),
                         curvature::FT             = FT(0.9),
                                qy::FT             = FT(0.4081632653061224)) where {FT}
    if r25==Inf
        r25 = model.VR.v_to_r * v25
    end
    r_leaf  = get_r(model.ReT, r25, t_leaf)
    vmax    = get_vmax(model.VcT, v25, t_leaf)
    jmax    = get_jmax(model.JT, j25, t_leaf)
    j       = get_j(jmax, par; curvature=curvature, qy=qy)
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
    get_an_ag_r_from_pi(; model, v25, p25, p_i, par, t_leaf, r25, qy)

Gross photosynthetic rate limited by carboxylation, given
- `model` A C4ParaSetVcJ type parameter sets that stores temperature denpendencies
- `v25` Maximal carboxylation rate at 298.15 K
- `p25` Maximal PEP carboxylation rate at 298.15 K
- `p_i` Leaf internal CO₂ partial pressure
- `par` Absorbed photosynthetic active radiation
- `t_leaf` Leaf temperature
- `r25` Respiration rate at 298.15 K
- `qy` Quantum yield of electron

The model is adapted from Collatz et al. (1992) "Coupled photosynthesis-stomatal conductance model for leaves of C4 plants."


"""
function get_an_ag_r_from_pi(;
                             model::C4ParaSetVcJ = C4VcVpJCLM{FT}(),
                               v25::FT           = FT(80.0),
                               p25::FT           = FT(100.0),
                               p_i::FT           = FT(30.0),
                               par::FT           = FT(1000.0),
                            t_leaf::FT           = FT(298.15),
                               r25::FT           = FT(Inf),
                                qy::FT           = FT(0.4081632653061224)) where {FT}
    if r25==Inf
        r25 = model.VR.v_to_r * v25
    end
    r_leaf  = get_r(model.ReT, r25, t_leaf)
    pmax    = get_vmax(model.VpT, p25, t_leaf)
    a_v     = get_vmax(model.VcT, v25, t_leaf)
    kpep    = get_kpep(model.kpT, t_leaf)
    # different a_j calculation from the C3 model
    a_j     = qy * par / 6
    a_p     = pmax * p_i / (p_i + kpep)
    ag      = min(a_v, a_j, a_p)
    an      = ag - r_leaf
    return an,ag,r_leaf
end
