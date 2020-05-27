"""
    get_empirical_gsw_from_model(model::ESMBallBerry, an::FT, p_atm::FT, p_i::FT, rh::FT, vpd::FT, Γ_star::FT)

Steady stage gsw from empirical approach given
- `model` An ESMBallBerry empirical model type parameter set
- `an` Net photosynthetic rate `[μmol m⁻² s⁻¹]`
- `p_atm` Atmospheric pressure
- `p_i` Leaf internal CO₂ partial pressure
- `rh` Relative humidity
- `vpd` Vapor pressure deficit in the air
- `Γ_star` CO₂ compensation point with the absence of dark respiration

Some parameters are not used for the Ball-Berry type model.

"""
function get_empirical_gsw_from_model(model::ESMBallBerry, an::FT, p_atm::FT, p_i::FT, rh::FT, vpd::FT, Γ_star::FT) where {FT}
    return model.g0 + model.g1 * rh * an / (p_i / p_atm * FT(1e6))
end




"""
    get_empirical_gsw_from_model(model::ESMGentine, an::FT, k_ratio::FT, p_atm::FT, p_i::FT)

Steady stage gsw from empirical approach given
- `model` An ESMGentine empirical model type parameter set
- `an` Net photosynthetic rate `[μmol m⁻² s⁻¹]`
- `p_atm` Atmospheric pressure
- `p_i` Leaf internal CO₂ partial pressure
- `rh` Relative humidity
- `vpd` Vapor pressure deficit in the air
- `Γ_star` CO₂ compensation point with the absence of dark respiration

Some parameters are not used for the Leuning type model.

"""
function get_empirical_gsw_from_model(model::ESMGentine, an::FT, k_ratio::FT, p_atm::FT, p_i::FT) where {FT}
    return model.g0 + model.g1 * an / (p_i / p_atm * FT(1e6)) * k_ratio
end




"""
    get_empirical_gsw_from_model(model::ESMLeuning, an::FT, p_atm::FT, p_i::FT, rh::FT, vpd::FT, Γ_star::FT)

Steady stage gsw from empirical approach given
- `model` An ESMLeuning empirical model type parameter set
- `an` Net photosynthetic rate `[μmol m⁻² s⁻¹]`
- `p_atm` Atmospheric pressure
- `p_i` Leaf internal CO₂ partial pressure
- `rh` Relative humidity
- `vpd` Vapor pressure deficit in the air
- `Γ_star` CO₂ compensation point with the absence of dark respiration

Some parameters are not used for the Leuning type model.

"""
function get_empirical_gsw_from_model(model::ESMLeuning, an::FT, p_atm::FT, p_i::FT, rh::FT, vpd::FT, Γ_star::FT) where {FT}
    return model.g0 + model.g1 * an / ((p_i - Γ_star) / p_atm * FT(1e6)) / (1 + vpd/model.d0)
end




"""
    get_empirical_gsw_from_model(model::ESMMedlyn, an::FT, p_atm::FT, p_i::FT, rh::FT, vpd::FT, Γ_star::FT)

Steady stage gsw from empirical approach given
- `model` An ESMMedlyn empirical model type parameter set
- `an` Net photosynthetic rate `[μmol m⁻² s⁻¹]`
- `p_atm` Atmospheric pressure
- `p_i` Leaf internal CO₂ partial pressure
- `rh` Relative humidity
- `vpd` Vapor pressure deficit in the air
- `Γ_star` CO₂ compensation point with the absence of dark respiration

Some parameters are not used for the Medlyn type model.

"""
function get_empirical_gsw_from_model(model::ESMMedlyn, an::FT, p_atm::FT, p_i::FT, rh::FT, vpd::FT, Γ_star::FT) where {FT}
    return model.g0 + (1 + model.g1/sqrt(vpd)) * an / (p_i / p_atm * FT(1e6))
end




"""
    get_empirical_gsw(curvature::FT, g_ias_c::FT, g_ias_e::FT, g_max::FT, j25::FT, p25::FT, p_a::FT, p_atm::FT, p_H₂O::FT, p_O₂::FT, par::FT, qy::FT, r25::FT, t_leaf::FT, v25::FT, paraset::C3ParaSet, model::AbstractEmpiricalStomatalModel)

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
- `paraset` An C3ParaSet type parameter set
- `model` An AbstractEmpiricalStomatalModel empirical model type


"""
function get_empirical_gsw(
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
    gsw_min = model.g0
    gsw_max = g_max * (t_leaf / FT(K_25))^FT(1.8)
    gsw     = (gsw_min + gsw_max) / 2
    rh      = p_H₂O / get_saturated_vapor_pressure(t_leaf)
    vpd     = get_saturated_vapor_pressure(t_leaf) - p_H₂O
    Γ_star  = get_Γ_star(paraset.ΓsT, t_leaf)
    while true
        gsw = (gsw_min + gsw_max) / 2
        gsc = gsw / FT(1.6) / ( 1 + g_ias_c * gsw^(g_ias_e) )
        an,ag,r,p_i = get_an_ag_r_pi_from_gsc(
                                              paraset,
                                              gsc,
                                              v25,
                                              j25,
                                              p25,
                                              p_a,
                                              t_leaf,
                                              par,
                                              p_atm,
                                              p_O₂,
                                              r25,
                                              curvature,
                                              qy)
        gsm = get_empirical_gsw_from_model(model, an, p_atm, p_i, rh, vpd, Γ_star)
        if (abs(gsm-gsw) < 1e-5) || (gsw_max-gsw_min<1e-5)
            break
        elseif gsm<gsw
            gsw_max = gsw
        else
            gsw_min = gsw
        end
    end

    # return the result
    return gsw
end




"""
    get_empirical_gsw(curvature::FT, g_ias_c::FT, g_ias_e::FT, g_max::FT, j25::FT, p25::FT, p_a::FT, p_atm::FT, p_H₂O::FT, p_O₂::FT, par::FT, qy::FT, r25::FT, t_leaf::FT, v25::FT, paraset::C3ParaSet, model::AbstractEmpiricalStomatalModel)

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
- `paraset` An C4ParaSet type parameter set
- `model` An AbstractEmpiricalStomatalModel empirical model type


"""
function get_empirical_gsw(
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
    gsw_min = model.g0
    gsw_max = g_max * (t_leaf / FT(K_25))^FT(1.8)
    gsw     = (gsw_min + gsw_max) / 2
    rh      = p_H₂O / get_saturated_vapor_pressure(t_leaf)
    vpd     = get_saturated_vapor_pressure(t_leaf) - p_H₂O
    Γ_star  = FT(0.0)
    while true
        gsw = (gsw_min + gsw_max) / 2
        gsc = gsw / FT(1.6) / ( 1 + g_ias_c * gsw^(g_ias_e) )
        an,ag,r,p_i = get_an_ag_r_pi_from_gsc(
                                              paraset,
                                              gsc,
                                              v25,
                                              j25,
                                              p25,
                                              p_a,
                                              t_leaf,
                                              par,
                                              p_atm,
                                              p_O₂,
                                              r25,
                                              curvature,
                                              qy)
        gsm = get_empirical_gsw_from_model(model, an, p_atm, p_i, rh, vpd, Γ_star)
        if (abs(gsm-gsw) < 1e-5) || (gsw_max-gsw_min<1e-5)
            break
        elseif gsm<gsw
            gsw_max = gsw
        else
            gsw_min = gsw
        end
    end

    # return the result
    return gsw
end
