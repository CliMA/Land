"""
    get_an_ag_r_pi_from_gsc(model::C3ParaSet, gsc::FT, v25::FT, j25::FT, p25::FT, p_a::FT, t_leaf::FT, par::FT, p_atm::FT, p_O₂::FT, r25::FT, curvature::FT, qy::FT)

Net photosynthetic rate `tar_an`, gross photosynthetic rate `tar_ag`, respiration rate `tar_r`, and leaf internal CO₂ partial pressure `tar_p`, given
- `model` A C3ParaSet type parameter set
- `gsc` Leaf diffusive conductance to CO₂
- `v25` Maixmal carboxylation rate at 298.15 K (25 Celcius)
- `j25` Maximal electron transport rate at 298.15 K (25 Celcius)
- `p25` Maximal PEP carboxylation rate at 298.15 K, C4 plant only
- `p_a` Atmospheric CO₂ partial pressure
- `t_leaf` Leaf temperature
- `par` Photosynthetic active radiation
- `p_atm` Atmospheric pressure
- `p_O₂` O₂ partial pressure
- `r25` Leaf respiration rate at 298.15 K (25 Celcius). If ==Inf, r25 will be computed from v25
- `curvature` Curvature parameter to calculate J
- `qy` Quantum yield of electron

For C3 plants, w/wo ATP limitation. Some parameters are useless in C3 photosynthesis, like p25.

"""
function get_an_ag_r_pi_from_gsc(model::C3ParaSet = C3VcJBernacchi{FT}(),
                                   gsc::FT        = FT(0.1),
                                   v25::FT        = FT(80.0),
                                   j25::FT        = FT(135.0),
                                   p25::FT        = FT(120.0),
                                   p_a::FT        = FT(40.0),
                                t_leaf::FT        = FT(298.15),
                                   par::FT        = FT(1000.0),
                                 p_atm::FT        = FT(101325.0),
                                  p_O₂::FT        = FT(21278.25),
                                   r25::FT        = FT(Inf),
                             curvature::FT        = FT(0.9),
                                    qy::FT        = FT(0.4081632653061224)) where {FT}
    # compute a_net using bi-section method
    tar_p  = FT(0.0)
    tar_ag = FT(0.0)
    tar_an = FT(0.0)
    tar_r  = FT(0.0)
    max_p  = p_a
    min_p  = FT(0.1)
    while true
        tar_p = FT(0.5) * (max_p+min_p)
        an,ag,r = get_an_ag_r_from_pi(model,
                                      tar_p,
                                      v25,
                                      j25,
                                      par,
                                      p_O₂,
                                      t_leaf,
                                      r25,
                                      curvature,
                                      qy)
        tmp_g = an * FT(1e-6) / (p_a-tar_p) * p_atm

        # increase min_p when g is smaller than target
        if abs(tmp_g-gsc) < 1e-5
            tar_ag = ag
            tar_an = an
            tar_r  = r
            break
        elseif tmp_g<gsc
            min_p = tar_p
        else
            max_p = tar_p
        end

        # if the max_p and min_p equals, break; used for the case when g=0
        if abs(max_p-min_p) < 1e-5
            tar_ag = ag
            tar_an = an
            tar_r  = r
            break
        end
    end

    # return the an, ag, r, and p_i
    return tar_an, tar_ag, tar_r, tar_p
end




"""
    get_an_ag_r_pi_from_gsc(model::C4ParaSet, gsc::FT, v25::FT, j25::FT, p25::FT, p_a::FT, t_leaf::FT, par::FT, p_atm::FT, p_O₂::FT, r25::FT, curvature::FT, qy::FT)

Net photosynthetic rate `tar_an`, gross photosynthetic rate `tar_ag`, respiration rate `tar_r`, and leaf internal CO₂ partial pressure `tar_p`, given
- `model` A C4ParaSet type parameter set
- `gsc` Leaf diffusive conductance to CO₂
- `v25` Maixmal carboxylation rate at 298.15 K (25 Celcius)
- `j25` Maximal electron transport rate at 298.15 K (25 Celcius)
- `p25` Maximal PEP carboxylation rate at 298.15 K, C4 plant only
- `p_a` Atmospheric CO₂ partial pressure
- `t_leaf` Leaf temperature
- `par` Photosynthetic active radiation
- `p_atm` Atmospheric pressure
- `p_O₂` O₂ partial pressure
- `r25` Leaf respiration rate at 298.15 K (25 Celcius). If ==Inf, r25 will be computed from v25
- `curvature` Curvature parameter to calculate J
- `qy` Quantum yield of electron

For C4 plants. Some parameters are useless for C4 photosynthesis, like j25, p_O₂, and curvature.

"""
function get_an_ag_r_pi_from_gsc(model::C4ParaSet = C4VcVpJCLM{FT}(),
                                   gsc::FT        = FT(0.1),
                                   v25::FT        = FT(80.0),
                                   j25::FT        = FT(135.0),
                                   p25::FT        = FT(120.0),
                                   p_a::FT        = FT(40.0),
                                t_leaf::FT        = FT(298.15),
                                   par::FT        = FT(1000.0),
                                 p_atm::FT        = FT(101325.0),
                                  p_O₂::FT        = FT(21278.25),
                                   r25::FT        = FT(Inf),
                             curvature::FT        = FT(0.9),
                                    qy::FT        = FT(0.4081632653061224)) where {FT}
    # compute a_net using bi-section method
    tar_p  = FT(0.0)
    tar_ag = FT(0.0)
    tar_an = FT(0.0)
    tar_r  = FT(0.0)
    max_p  = p_a
    min_p  = FT(0.1)
    while true
        tar_p = FT(0.5) * (max_p+min_p)
        an,ag,r = get_an_ag_r_from_pi(model,
                                      tar_p,
                                      v25,
                                      p25,
                                      par,
                                      t_leaf,
                                      r25,
                                      qy)
        tmp_g = an * FT(1e-6) / (p_a-tar_p) * p_atm

        # increase min_p when g is smaller than target
        if abs(tmp_g-gsc) < 1e-5
            tar_ag = ag
            tar_an = an
            tar_r  = r
            break
        elseif tmp_g<gsc
            min_p = tar_p
        else
            max_p = tar_p
        end

        # if the max_p and min_p equals, break; used for the case when g=0
        if abs(max_p-min_p) < 1e-5
            tar_ag = ag
            tar_an = an
            tar_r  = r
            break
        end
    end

    # return the an, ag, r, and p_i
    return tar_an, tar_ag, tar_r, tar_p
end
