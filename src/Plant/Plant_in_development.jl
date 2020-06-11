###############################################################################
#
# Calculate Weibull function considering temperature
#
###############################################################################
function weibull_k_ratio(b::FT, c::FT, p::FT, T::FT, p_his::FT) where {FT}
    p_25  = min( p_his, p / relative_surface_tension(T) )
    kr_25 = exp( -1 * (-p_25/b) ^ (c) ) / relative_viscosity(T)
    return kr_25
end

function weibull_k_ratio(b::Array, c::Array, p::Array, T::Array, p_his::Array)
    FT = eltype(T)

    p_25  = min.( p_his, p ./ relative_surface_tension(T) )
    kr_25 = exp.( -1 .* (-p_25 ./ b) .^ c ) ./ relative_viscosity(T)
    return kr_25
end

function weibull_k_leaf_ratio(leaf::Leaf{FT}, T::FT) where {FT}
    p_25  = leaf.p_element[end] / relative_surface_tension(T)
    kr_25 = exp( -1 * (-p_25 / leaf.b) ^ (leaf.c) ) / relative_viscosity(T)
    return kr_25
end

function weibull_k_leaf_ratio(leaves::Array, T::Array)
    FT     = eltype(T)
    p_list = [leaves[i].p_element[end] for i in 1:length(T)]
    b_list = [leaves[i].b              for i in 1:length(T)]
    c_list = [leaves[i].c              for i in 1:length(T)]
    p_25   = p_list ./ relative_surface_tension(T)
    kr_25  = exp.( -1 .* (-p_list ./ b_list) .^ c_list ) ./ relative_viscosity(T)
    return kr_25
end








# Move to Leaf?
"""
    get_empirical_gsw_pi(curvature::FT, g_ias_c::FT, g_ias_e::FT, g_max::FT, j25::FT, p25::FT, p_a::FT, p_atm::FT, p_H₂O::FT, p_O₂::FT, par::FT, qy::FT, r25::FT, t_leaf::FT, v25::FT, paraset::C3ParaSet, model::AbstractEmpiricalStomatalModel)
    get_empirical_gsw_pi(curvature::FT, g_ias_c::FT, g_ias_e::FT, g_max::FT, j25::FT, p25::FT, p_a::FT, p_atm::FT, p_H₂O::FT, p_O₂::FT, par::FT, qy::FT, r25::FT, t_leaf::FT, v25::FT, paraset::C3ParaSet, model::AbstractEmpiricalStomatalModel)

Steady state gsw from empirical approach, given
- `curvature` Curvature parameter to calculate J
- `g_ias_c` Mesophyll conductance correction multiplier
- `g_ias_e` Mesophyll conductance correction exponent
- `g_max` Maximal leaf diffusive conductance at 298.15 K
- `j25` Maximal electron transport rate at 298.15 K (25 Celcius)
- `p25` Maximal PEP carboxylation rate at 298.15 K, C4 plant only
- `p_a` Atmospheric CO₂ partial pressure
- `p_atm` Atmospheric pressure
- `p_H₂O` Atmospheric H₂O partial pressure
- `p_O₂` O₂ partial pressure
- `par` Photosynthetic active radiation
- `qy` Quantum yield of electron
- `r25` Leaf respiration rate at 298.15 K (25 Celcius). If ==Inf, r25 will be computed from v25
- `t_leaf` Leaf temperature
- `v25` Maixmal carboxylation rate at 298.15 K (25 Celcius)
- `paraset` An C3ParaSet or C4paraSet type parameter set
- `model` An AbstractEmpiricalStomatalModel empirical model type
"""
function get_empirical_gsw_pi(
                         curvature::FT,
                           g_ias_c::FT,
                           g_ias_e::FT,
                             g_max::FT,
                               j25::FT,
                               p25::FT,
                               p_a::FT,
                             p_atm::FT,
                             p_H₂O::FT,
                              p_O₂::FT,
                               par::FT,
                                qy::FT,
                               r25::FT,
                            t_leaf::FT,
                               v25::FT,
                           paraset::C3ParaSet,
                             model::AbstractEmpiricalStomatalModel) where {FT}
    gsw_min = model.g0 * (t_leaf / FT(K_25))^FT(1.8)
    gsw_max = g_max * (t_leaf / FT(K_25))^FT(1.8)
    rh      = p_H₂O / saturation_vapor_pressure(t_leaf)
    vpd     = saturation_vapor_pressure(t_leaf) - p_H₂O
    Γ_star  = get_Γ_star(paraset.ΓsT, t_leaf)

    # Use iteration to obtain target p_i
    gsw   = 0
    p_i   = p_a - 1
    count = 0
    while true
        count += 1
        if count > 50
            printstyled("Total iteration exceeds 50 times for get_empirical_gsw function.\n", color=:red)
            break
        end

        an,ag,r = an_ag_r_from_pi(paraset,
                                  p_i,
                                  v25,
                                  j25,
                                  par,
                                  p_O₂,
                                  t_leaf,
                                  r25,
                                  curvature,
                                  qy)
        g_aci_0  = an*FT(1e-6) / (p_a - p_i) * p_atm
        g_aci_0 *= FT(1.6) * ( 1 + g_ias_c * gsw^(g_ias_e) )
        g_mod_0  = get_empirical_gsw_from_model(model, an, p_atm, p_i, rh, vpd, Γ_star, FT(1.0))
        g_diff_0 = g_aci_0 - g_mod_0

        if g_diff_0 < 1e-4
            gsw = max(gsw_min, min(g_aci_0, gsw_max))
            break
        end

        p_j = p_i + FT(0.01)
        ao,ag,r = an_ag_r_from_pi(paraset,
                                  p_i,
                                  v25,
                                  j25,
                                  par,
                                  p_O₂,
                                  t_leaf,
                                  r25,
                                  curvature,
                                  qy)
        g_aci_1  = ao*FT(1e-6) / (p_a - p_j) * p_atm
        g_aci_1 *= FT(1.6) * ( 1 + g_ias_c * gsw^(g_ias_e) )
        g_mod_1  = get_empirical_gsw_from_model(model, ao, p_atm, p_j, rh, vpd, Γ_star, FT(1.0))
        g_diff_1 = g_aci_1 - g_mod_1

        slope = (g_diff_1 - g_diff_0) * 100
        p_i  -= g_diff_0 / slope
    end

    # return the result
    return gsw,p_i
end

function get_empirical_gsw_pi(
                         curvature::FT,
                           g_ias_c::FT,
                           g_ias_e::FT,
                             g_max::FT,
                               j25::FT,
                               p25::FT,
                               p_a::FT,
                             p_atm::FT,
                             p_H₂O::FT,
                              p_O₂::FT,
                               par::FT,
                                qy::FT,
                               r25::FT,
                            t_leaf::FT,
                               v25::FT,
                           paraset::C4ParaSet,
                             model::AbstractEmpiricalStomatalModel) where {FT}
    gsw_min = model.g0 * (t_leaf / FT(K_25))^FT(1.8)
    gsw_max = g_max * (t_leaf / FT(K_25))^FT(1.8)
    rh      = p_H₂O / saturation_vapor_pressure(t_leaf)
    vpd     = saturation_vapor_pressure(t_leaf) - p_H₂O
    Γ_star  = FT(0.0)

    # Use iteration to obtain target p_i
    gsw   = 0
    p_i   = p_a - 1
    count = 0
    while true
        count += 1
        if count > 50
            printstyled("Total iteration exceeds 50 times for get_empirical_gsw function.\n", color=:red)
            break
        end

        an,ag,r = an_ag_r_from_pi(paraset,
                                  p_i,
                                  v25,
                                  p25,
                                  par,
                                  t_leaf,
                                  r25,
                                  qy)
        g_aci_0  = an*FT(1e-6) / (p_a - p_i) * p_atm
        g_aci_0 *= FT(1.6) * ( 1 + g_ias_c * gsw^(g_ias_e) )
        g_mod_0  = get_empirical_gsw_from_model(model, an, p_atm, p_i, rh, vpd, Γ_star, FT(1.0))
        g_diff_0 = g_aci_0 - g_mod_0

        if g_diff_0 < 1e-4
            gsw = max(gsw_min, min(g_aci_0, gsw_max))
            break
        end

        p_j = p_i + FT(0.01)
        ao,ag,r = an_ag_r_from_pi(paraset,
                                  p_j,
                                  v25,
                                  p25,
                                  par,
                                  t_leaf,
                                  r25,
                                  qy)
        g_aci_1  = ao*FT(1e-6) / (p_a - p_j) * p_atm
        g_aci_1 *= FT(1.6) * ( 1 + g_ias_c * gsw^(g_ias_e) )
        g_mod_1  = get_empirical_gsw_from_model(model, ao, p_atm, p_j, rh, vpd, Γ_star, FT(1.0))
        g_diff_1 = g_aci_1 - g_mod_1

        slope = (g_diff_1 - g_diff_0) * 100
        p_i  -= g_diff_0 / slope
    end

    # return the result
    return gsw,p_i
end




















# TODO Move to Phtosynthesis module
"""
    get_a_par_curve(v25, j25, Γ_star, gsc, p_a, tem, p_O₂, r25)

Lists of phtosynthetic active radiation `list_par`, net phtosynthetic rate `list_an`, gross phtosynthetic rate `list_ag`, and leaf internal CO₂ partial pressure `list_pi`, given constant
- `v25` Maximal carboxylation rate at 298.15 K (25 Celcius)
- `j25` Maximal electron transport rate at 298.15 K
- `Γ_star` CO₂ compensation point with the absence of dark respiration
- `gsc` Leaf diffusive conductance to CO₂
- `p_a` Atmospheric CO₂ partial pressure
- `tem` Leaf temperature
- `p_atm` Atmospheric pressure
- `p_O₂` O₂ partial pressure
- `r25` Leaf respiration rate at 298.15 K. If ==Inf, use v25 to compute r25.

# Description
This function is meant to obtain the Anet-PAR curve.
This function may also be used to test the phtosynthesis module.
"""
function get_a_par_curve(model::AbstractPhotoModelParaSet,
                           v25::FT,
                           j25::FT,
                           p25::FT,
                           gsc::FT,
                           p_a::FT,
                           tem::FT,
                         p_atm::FT,
                          p_O₂::FT,
                           r25::FT,
                     curvature::FT,
                            qy::FT) where {FT}
    # generate a list of p_i from gamma to 200 Pa
    list_par = [FT(val) for val in 5.0:5.0:2000.0]
    list_ag  = list_par .* FT(-Inf)
    list_an  = list_par .* FT(-Inf)
    list_pi  = list_par .* FT(-Inf)

    # iterate through the list
    for i in 1:length(list_par)
        par         = list_par[i]
        an,ag,r,p_i = an_ag_r_pi_from_gsc(
                                          model,
                                          gsc,
                                          v25,
                                          j25,
                                          p25,
                                          p_a,
                                          tem,
                                          par,
                                          p_atm,
                                          p_O₂,
                                          r25,
                                          curvature,
                                          qy)
        list_ag[i]  = ag
        list_an[i]  = an
        list_pi[i]  = p_i
    end

    # return the arrays
    return list_par,list_an,list_ag,list_pi
end






"""
    get_a_pi_curve(v25, j25, Γ_star, tem, par, p_O₂, r25)

Lists of leaf internal CO₂ partial pressure `list_pi`, net photosynthetic rate `list_an`, and gross photosynthetic rate `list_ag`, given constant
- `v25` Maximal carboxylation rate at 298.15 K (25 Celcius)
- `j25` Maximal electron transport rate at 298.15 K
- `Γ_star` CO₂ compensation point with the absence of dark respiration
- `tem` Leaf temperature
- `par` Photosynthetic active radiation
- `p_O₂` O₂ partial pressure
- `r25` Leaf respiration rate at 298.15 K. If ==Inf, use v25 to compute r25.
"""
function get_a_pi_curve(model::C3ParaSet,
                          v25::FT,
                          j25::FT,
                          p25::FT,
                          tem::FT,
                          par::FT,
                         p_O₂::FT,
                          r25::FT,
                    curvature::FT,
                           qy::FT) where {FT}
    # generate a list of p_i from gamma to 200 Pa
    list_pi = [FT(val) for val in 1.0:1.0:200.0]
    list_ag = list_pi .* FT(-Inf)
    list_an = list_pi .* FT(-Inf)

    # iterate through the list
    for i in 1:length(list_pi)
        p_i        = list_pi[i]
        a_n,a_g,r  = an_ag_r_from_pi(
                                     model,
                                     p_i,
                                     v25,
                                     j25,
                                     par,
                                     p_O₂,
                                     tem,
                                     r25,
                                     curvature,
                                     qy)
        list_ag[i] = a_g
        list_an[i] = a_n
    end

    # return the arrays
    return list_pi,list_an,list_ag
end

