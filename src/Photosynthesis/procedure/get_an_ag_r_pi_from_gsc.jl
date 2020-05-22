"""
    get_an_ag_r_pi_from_gsc(; model, v25, j25, gsc, p_a, t_leaf, par, p_atm, p_O₂, r25, curvature, qy)

Net photosynthetic rate `tar_an`, gross photosynthetic rate `tar_ag`, respiration rate `tar_r`, and leaf internal CO₂ partial pressure `tar_p`, given
- `model` A C3ParaSet type parameter set
- `v25` Maixmal carboxylation rate at 298.15 K (25 Celcius)
- `j25` Maximal electron transport rate at 298.15 K (25 Celcius)
- `gsc` Leaf diffusive conductance to CO₂
- `p_a` Atmospheric CO₂ partial pressure
- `t_leaf` Leaf temperature
- `par` Photosynthetic active radiation
- `p_atm` Atmospheric pressure
- `p_O₂` O₂ partial pressure
- `r25` Leaf respiration rate at 298.15 K (25 Celcius). If ==Inf, r25 will be computed from v25
- `curvature` Curvature parameter to calculate J
- `qy` Quantum yield of electron

"""
function get_an_ag_r_pi_from_gsc(;
                                 model::C3ParaSet = C3VcJBernacchi{FT}(),
                                   v25::FT        = FT(80.0),
                                   j25::FT        = FT(135.0),
                                   gsc::FT        = FT(0.1),
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
        an,ag,r = get_an_ag_r_from_pi(
                                      model = model,
                                        v25 = v25,
                                        j25 = j25,
                                        p_i = tar_p,
                                        par = par,
                                       p_O₂ = p_O₂,
                                     t_leaf = t_leaf,
                                        r25 = r25,
                                  curvature = curvature,
                                         qy = qy)
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
    get_an_ag_r_pi_from_gsc(; model, v25, p25, gsc, p_a, t_leaf, par, p_atm, r25, qy)

Net photosynthetic rate `tar_an`, gross photosynthetic rate `tar_ag`, respiration rate `tar_r`, and leaf internal CO₂ partial pressure `tar_p`, given
- `model` A C4ParaSet type parameter set
- `v25` Maixmal carboxylation rate at 298.15 K (25 Celcius)
- `p25` Maximal PEP carboxylation rate at 298.15 K (25 Celcius)
- `gsc` Leaf diffusive conductance to CO₂
- `p_a` Atmospheric CO₂ partial pressure
- `t_leaf` Leaf temperature
- `par` Photosynthetic active radiation
- `p_atm` Atmospheric pressure
- `r25` Leaf respiration rate at 298.15 K (25 Celcius). If ==Inf, r25 will be computed from v25
- `qy` Quantum yield of electron


"""
function get_an_ag_r_pi_from_gsc(;
                                 model::C4ParaSet = C4VcVpJCLM{FT}(),
                                   v25::FT        = FT(80.0),
                                   p25::FT        = FT(120.0),
                                   gsc::FT        = FT(0.1),
                                   p_a::FT        = FT(40.0),
                                t_leaf::FT        = FT(298.15),
                                   par::FT        = FT(1000.0),
                                 p_atm::FT        = FT(101325.0),
                                   r25::FT        = FT(Inf),
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
        an,ag,r = get_an_ag_r_from_pi(
                                      model = model,
                                        v25 = v25,
                                        p25 = p25,
                                        p_i = tar_p,
                                        par = par,
                                     t_leaf = t_leaf,
                                        r25 = r25,
                                         qy = qy)
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
