###############################################################################
#
# Arrhenius corrections
# These functions passed the FT test
# These functions are documented in the Photosynthesis page
#
###############################################################################
"""
    arrhenius_correction(td_set::ArrheniusTD, T::FT)
    arrhenius_correction(td_set::ArrheniusPeakTD, T::FT)

A correction factor based on arrhenius's fitting procedure, given
- `td_set` An `ArrheniusTD` or `ArrheniusPeakTD` type struct
- `T` Leaf temperature in `[K]`

The equation used for `ArrheniusTD` is `correction = exp( c - ΔHa/(R*T_leaf) )`.

The equations used for `ArrheniusPeakTD` are
```
    C = 1 + exp[(Sv*T_0 - Hd )/(RT_0 )]
    f_above = C * exp [ ( Ha /RT_0 ) * ( 1 - T_0/T_1 ) ]
    f_below = 1 + exp [ ( Sv*T_1 - Hd ) / ( RT_1) ]
    correction = f_above /f_below
```
"""
function arrhenius_correction(td_set::ArrheniusTD, T::FT) where {FT}
    return exp( -td_set.ΔHa_to_R/T + td_set.ΔHa_to_RT25 )
end

function arrhenius_correction(td_set::ArrheniusTD, T::Array)
    return exp.( -td_set.ΔHa_to_R ./ T .+ td_set.ΔHa_to_RT25 )
end

function arrhenius_correction(td_set::ArrheniusPeakTD, T::FT) where {FT}
    @unpack C, ΔHa_to_RT25, ΔHd_to_R, ΔSv_to_R = td_set
    return C * exp( ΔHa_to_RT25 * (1-FT(K_25)/T) ) / ( 1 + exp(ΔSv_to_R - ΔHd_to_R/T) )
end

function arrhenius_correction(td_set::ArrheniusPeakTD, T::Array)
    @unpack C, ΔHa_to_RT25, ΔHd_to_R, ΔSv_to_R = td_set
    FT = eltype(T)
    return C .* exp.( ΔHa_to_RT25 .* (1 .- FT(K_25) ./ T) ) ./ ( 1 .+ exp.(ΔSv_to_R .- ΔHd_to_R ./ T) )
end








###############################################################################
#
# Temperature dependency of the photosynthetic parameters
# These functions passed the FT test
# These functions are documented in the Photosynthesis page
#
###############################################################################
"""
    get_jmax(td_set::AbstractJmaxTD, j25::FT, T::FT)

Maximal electron transport rate at leaf temperature, given
- `td_set` AbstractJmax parameter set
- `j25` Maximal electron transport rate at 298.15 K
- `T` Leaf temperature
"""
function get_jmax(td_set::AbstractTDParameterSet, j25::FT, T::FT) where {FT}
    return j25 * arrhenius_correction(td_set, T)
end




"""
    get_kc(td_set::AbstractKc, T::FT)

Kc at leaf temperature, given
- `td_set` One `AbstractKc` type that store temperature correction information
- `T` Leaf temperature
"""
function get_kc(td_set::ArrheniusTD, T::FT) where {FT}
    return td_set.VAL_25 * arrhenius_correction(td_set, T)
end




"""
    get_ko(td_set::AbstractKo, T::FT)

Ko at leaf temperature, given
- `td_set` One `AbstractKo` type that store temperature correction information
- `T` Leaf temperature
"""
function get_ko(td_set::ArrheniusTD, T::FT) where {FT}
    return td_set.VAL_25 * arrhenius_correction(td_set, T)
end




"""
    get_kpep(td_set::AbstractKpep, T::FT)

Kpep at leaf temperature, given
- `td_set` One `AbstractKpep` type that store temperature correction information
- `T` Leaf temperature
"""
function get_kpep(td_set::ArrheniusTD, T::FT) where {FT}
    return td_set.VAL_25 * arrhenius_correction(td_set, T)
end




"""
    get_r(td_set::AbstractRespirationTD, r25::FT, T::FT)

Leaf respiration rate `r_leaf` at leaf temperature, given
- `td_set` One `AbstractRespirationTD` type that store temperature correction information
- `r25` Leaf respiration rate at 298.15 K (25 Celcius)
- `T` Leaf temperature
"""
function get_r(td_set::AbstractTDParameterSet, r25::FT, T::FT) where {FT}
    return r25 * arrhenius_correction(td_set, T)
end




"""
    get_vmax(v25::FT, td_set::AbstractVmaxTD, T::FT)

Maximal electron transport rate at leaf temperature, given
- `td_set` AbstractVmaxTD parameter set
- `v25` Maximal carboxylation rate at 298.15 K
- `T` Leaf temperature
"""
function get_vmax(td_set::AbstractTDParameterSet, v25::FT, T::FT) where {FT}
    return v25 * arrhenius_correction(td_set, T)
end




"""
    get_Γ_star(td_set::AbstractΓStarTD, T::FT)

Γ_star at leaf temperature, given
- `td_set` One `AbstractΓStarTD` type that store temperature correction information
- `T` Leaf temperature
"""
function get_Γ_star(td_set::ArrheniusTD, T::FT) where {FT}
    return td_set.VAL_25 * arrhenius_correction(td_set, T)
end

function get_Γ_star(td_set::ArrheniusTD, T::Array)
    return td_set.VAL_25 .* arrhenius_correction(td_set, T)
end








###############################################################################
#
# Calculate the electron transport rate
# These functions passed the FT test
# This function is documented in the Photosynthesis page
#
###############################################################################
"""
    get_j(jmax::FT, par::FT, curvature::FT, qy::FT)

Electron transport rate `j`, given
- `jmax` Maximal eclectron transport @ leaf temperature (not 298.15 K)
- `par` Absorbed PAR (photosynthetic active radiation)
- `curvature` Curvature parameter (default at 0.9)
- `qy` Quantum yield of electron, `qy = maxPSII * PSII_frac`, maxPSII is maximal PSII yield (default at 4.0/4.9), PSII_frac is Fraction of absorbed light used by PSII ETR (default at 0.5)
"""
function get_j(jmax::FT, par::FT, curvature::FT, qy::FT) where {FT}
    # a = curvature
    b = qy * par + jmax
    c = qy * par * jmax
    j = ( b - sqrt(b^2 - 4*curvature*c) ) / (2*curvature)
    return j
end








###############################################################################
#
# Calculate photosynthesis from CO₂ partial pressure
# These functions passed the FT test
# These functions are documented in the Photosynthesis page
#
###############################################################################
"""
    an_ag_r_from_pi(model::C3ParaSet,
                    p_i::FT,
                    v25::FT,
                    j25::FT,
                    par::FT,
                    p_O₂::FT,
                    T::FT,
                    r25::FT,
                    curvature::FT,
                    qy::FT)

Gross photosynthetic rate limited by carboxylation with ATP limitation, given
- `model` A C3ParaSet type parameter sets that stores temperature denpendencies
- `p_i` Leaf internal CO₂ partial pressure
- `v25` Maximal carboxylation rate at 298.15 K
- `j25` Maximal electron transport rate
- `par` Absorbed photosynthetic active radiation
- `p_O₂` O₂ partial pressure
- `T` Leaf temperature
- `r25` Respiration rate at 298.15 K
- `curvature` Curvature parameter to calculate J
- `qy` Quantum yield of electron

The equations used are from Farquhar et al. (1980) "A biochemical model of photosynthetic CO₂ assimilation in leaves of C3 species."

"""
function an_ag_r_from_pi(model::C3ParaSet,
                           p_i::FT,
                           v25::FT,
                           j25::FT,
                           par::FT,
                          p_O₂::FT,
                             T::FT,
                           r25::FT,
                     curvature::FT,
                            qy::FT) where {FT}
    if r25==Inf
        r25 = model.VR * v25
    end
    r_leaf  = get_r(model.ReT, r25, T)
    vmax    = get_vmax(model.VcT, v25, T)
    jmax    = get_jmax(model.JT, j25, T)
    j       = get_j(jmax, par, curvature, qy)
    kc      = get_kc(model.KcT, T)
    ko      = get_ko(model.KoT, T)
    km      = kc * (1 + p_O₂/ko)
    Γ_star  = get_Γ_star(model.ΓsT, T)
    a_p     = vmax / 2
    a_v     = vmax * (p_i-Γ_star) / (p_i+km)
    a_j     = j * (p_i-Γ_star) / ( 4*(p_i + 2*Γ_star) )
    ag      = min(a_p, a_v, a_j)
    an      = ag - r_leaf
    return an,ag,r_leaf
end




"""
    an_ag_r_from_pi(model::C4ParaSet,
                    p_i::FT,
                    v25::FT,
                    p25::FT,
                    par::FT,
                    T::FT,
                    r25::FT,
                    qy::FT)

Gross photosynthetic rate limited by carboxylation, given
- `model` A `C4ParaSet` type parameter sets that stores temperature denpendencies
- `p_i` Leaf internal CO₂ partial pressure
- `v25` Maximal carboxylation rate at 298.15 K
- `p25` Maximal PEP carboxylation rate at 298.15 K
- `par` Absorbed photosynthetic active radiation
- `T` Leaf temperature
- `r25` Respiration rate at 298.15 K
- `qy` Quantum yield of electron

The model is adapted from Collatz et al. (1992) "Coupled photosynthesis-stomatal conductance model for leaves of C4 plants."

"""
function an_ag_r_from_pi(model::C4ParaSet,
                           p_i::FT,
                           v25::FT,
                           p25::FT,
                           par::FT,
                             T::FT,
                           r25::FT,
                            qy::FT) where {FT}
    if r25==Inf
        r25 = model.VR * v25
    end
    r_leaf  = get_r(model.ReT, r25, T)
    pmax    = get_vmax(model.VpT, p25, T)
    a_v     = get_vmax(model.VcT, v25, T)
    kpep    = get_kpep(model.KpT, T)
    # different a_j calculation from the C3 model
    a_j     = qy * par / 6
    a_p     = pmax * p_i / (p_i + kpep)
    ag      = min(a_v, a_j, a_p)
    an      = ag - r_leaf
    return an,ag,r_leaf
end








###############################################################################
#
# Calculate photosynthesis from leaf diffusive conductance for CO₂
# These functions passed the FT test
# These functions are documented in the Photosynthesis page
#
###############################################################################
"""
    an_ag_r_pi_from_gsc(model::AbstractPhotoModelParaSet,
                        gsc::FT,
                        v25::FT,
                        j25::FT,
                        p25::FT,
                        p_a::FT,
                        T::FT,
                        par::FT,
                        p_atm::FT,
                        p_O₂::FT,
                        r25::FT,
                        curvature::FT,
                        qy::FT)

Net photosynthetic rate `tar_an`, gross photosynthetic rate `tar_ag`, respiration rate `tar_r`, and leaf internal CO₂ partial pressure `tar_p`, given
- `model` A `C3ParaSet` or `C4ParaSet` type parameter set
- `gsc` Leaf diffusive conductance to CO₂
- `v25` Maixmal carboxylation rate at 298.15 K (25 Celcius)
- `j25` Maximal electron transport rate at 298.15 K (25 Celcius)
- `p25` Maximal PEP carboxylation rate at 298.15 K, C4 plant only
- `p_a` Atmospheric CO₂ partial pressure
- `T` Leaf temperature
- `par` Photosynthetic active radiation
- `p_atm` Atmospheric pressure
- `p_O₂` O₂ partial pressure
- `r25` Leaf respiration rate at 298.15 K (25 Celcius). If ==Inf, r25 will be computed from v25
- `curvature` Curvature parameter to calculate J
- `qy` Quantum yield of electron

For C3 plants, with ATP limitation. Some parameters are useless in C3 photosynthesis, like p25.

For C4 plants. Some parameters are useless for C4 photosynthesis, like j25, p_O₂, and curvature.

"""
function an_ag_r_pi_from_gsc(model::C3ParaSet,
                               gsc::FT,
                               v25::FT,
                               j25::FT,
                               p25::FT,
                               p_a::FT,
                                 T::FT,
                               par::FT,
                             p_atm::FT,
                              p_O₂::FT,
                               r25::FT,
                         curvature::FT,
                                qy::FT) where {FT}
    # when there is not light
    if par < 1e-3
        if r25==Inf
            r25 = model.VR * v25
        end
        tar_r  = get_r(model.ReT, r25, T)
        tar_an = -tar_r
        tar_ag = FT(0.0)
        tar_pi = p_a + tar_r*FT(1e-6) * p_atm / gsc
    else
        int_y  = p_a / p_atm * gsc
        tar_pi = FT(0.1)
        tar_an = FT(0.0)
        tar_ag = FT(0.0)
        tar_r  = FT(0.0)
        count  = 0
        while true
            count  += 1
            an,ag,r = an_ag_r_from_pi(model,
                                      tar_pi,
                                      v25,
                                      j25,
                                      par,
                                      p_O₂,
                                      T,
                                      r25,
                                      curvature,
                                      qy)
            tar_g   = (int_y - an*FT(1e-6)) / tar_pi * p_atm

            if abs(tar_g - gsc) < 1e-6
                tar_an = an
                tar_ag = ag
                tar_r  = r
                break
            end

            if count > 50
                tar_an = an
                tar_ag = ag
                tar_r  = r
                printstyled("Total number of iterations exceeds 50 times for C3ParaSet with FT!\n", color=:red)
                println("\tThe provided PAR is ", par)
                println("\tThe provided gsc is ", gsc)
                println("\tThe provided V25 is ", v25)
                println("")
                break
            end

            tar_qi  = tar_pi + FT(0.001)
            ao,ag,r = an_ag_r_from_pi(model,
                                      tar_qi,
                                      v25,
                                      j25,
                                      par,
                                      p_O₂,
                                      T,
                                      r25,
                                      curvature,
                                      qy)
            tar_h   = (int_y - ao*FT(1e-6)) / tar_qi * p_atm

            slope   = (tar_h - tar_g) * 1000
            tar_pi += (gsc - tar_g) / slope
        end
    end

    # return the an, ag, r, and p_i
    return tar_an, tar_ag, tar_r, tar_pi
end

function an_ag_r_pi_from_gsc(model::C3ParaSet,
                               gsc::Array,
                               v25::FT,
                               j25::FT,
                               p25::FT,
                               p_a::FT,
                                 T::Array,
                               par::Array,
                             p_atm::FT,
                              p_O₂::FT,
                               r25::FT,
                         curvature::FT,
                                qy::FT) where {FT}
    # create lists to return
    N       = length(gsc)
    list_an = zeros(FT,N)
    list_ag = zeros(FT,N)
    list_r  = zeros(FT,N)
    list_pi = zeros(FT,N)

    # iterate through the list
    for i in 1:N
        _gsc = gsc[i]
        _par = par[i]
        _T   = T[i]

        # when there is no light
        if _par < 1e-3
            if r25==Inf
                r25 = model.VR * v25
            end
            tar_r  = get_r(model.ReT, r25, _T)
            tar_an = -tar_r
            tar_ag = FT(0.0)
            tar_pi = p_a + tar_r*FT(1e-6) * p_atm / _gsc
        else
            int_y  = p_a / p_atm * _gsc
            tar_pi = FT(0.1)
            tar_an = FT(0.0)
            tar_ag = FT(0.0)
            tar_r  = FT(0.0)
            count  = 0
            while true
                count  += 1
                an,ag,r = an_ag_r_from_pi(model,
                                          tar_pi,
                                          v25,
                                          j25,
                                          _par,
                                          p_O₂,
                                          _T,
                                          r25,
                                          curvature,
                                          qy)
                tar_g   = (int_y - an*FT(1e-6)) / tar_pi * p_atm

                if abs(tar_g - _gsc) < 1e-6
                    tar_an = an
                    tar_ag = ag
                    tar_r  = r
                    break
                end

                if count > 50
                    tar_an = an
                    tar_ag = ag
                    tar_r  = r
                    printstyled("Total number of iterations exceeds 50 times for C3ParaSet with Array!\n", color=:red)
                    println("\tThe provided PAR is ", par)
                    println("\tThe provided gsc is ", gsc)
                    println("\tThe provided V25 is ", v25)
                    println("")
                    break
                end

                tar_qi  = tar_pi + FT(0.001)
                ao,ag,r = an_ag_r_from_pi(model,
                                          tar_qi,
                                          v25,
                                          j25,
                                          _par,
                                          p_O₂,
                                          _T,
                                          r25,
                                          curvature,
                                          qy)
                tar_h   = (int_y - ao*FT(1e-6)) / tar_qi * p_atm

                slope   = (tar_h - tar_g) * 1000
                tar_pi += (_gsc - tar_g) / slope
            end
        end

        list_an[i] = tar_an
        list_ag[i] = tar_ag
        list_r[i]  = tar_r
        list_pi[i] = tar_pi
    end

    # return the an, ag, r, and p_i
    return list_an, list_ag, list_r, list_pi
end

function an_ag_r_pi_from_gsc(model::C4ParaSet,
                               gsc::FT,
                               v25::FT,
                               j25::FT,
                               p25::FT,
                               p_a::FT,
                                 T::FT,
                               par::FT,
                             p_atm::FT,
                              p_O₂::FT,
                               r25::FT,
                         curvature::FT,
                                qy::FT) where {FT}
    # when there is not light
    if par < 1e-3
        if r25==Inf
            r25 = model.VR * v25
        end
        tar_r  = get_r(model.ReT, r25, T)
        tar_an = -tar_r
        tar_ag = FT(0.0)
        tar_pi = p_a + tar_r*FT(1e-6) * p_atm / gsc
    else
        int_y  = p_a / p_atm * gsc
        tar_pi = FT(0.1)
        tar_an = FT(0.0)
        tar_ag = FT(0.0)
        tar_r  = FT(0.0)
        count  = 0
        while true
            count  += 1
            an,ag,r = an_ag_r_from_pi(model,
                                      tar_pi,
                                      v25,
                                      p25,
                                      par,
                                      T,
                                      r25,
                                      qy)
            tar_g   = (int_y - an*FT(1e-6)) / tar_pi * p_atm

            if abs(tar_g - gsc) < 1e-6
                tar_an = an
                tar_ag = ag
                tar_r  = r
                break
            end

            if count > 50
                tar_an = an
                tar_ag = ag
                tar_r  = r
                printstyled("Total number of iterations exceeds 50 times for C4ParaSet with FT!\n", color=:red)
                println("\tThe provided PAR is ", par)
                println("\tThe provided gsc is ", gsc)
                println("\tThe provided V25 is ", v25)
                println("")
                break
            end

            tar_q   = tar_pi + FT(0.001)
            ao,ag,r = an_ag_r_from_pi(model,
                                      tar_q,
                                      v25,
                                      p25,
                                      par,
                                      T,
                                      r25,
                                      qy)
            tar_h   = (int_y - ao*FT(1e-6)) / tar_q * p_atm

            slope   = (tar_h - tar_g) * 1000
            tar_pi += (gsc - tar_g) / slope
        end
    end

    # return the an, ag, r, and p_i
    return tar_an, tar_ag, tar_r, tar_pi
end

function an_ag_r_pi_from_gsc(model::C4ParaSet,
                               gsc::Array,
                               v25::FT,
                               j25::FT,
                               p25::FT,
                               p_a::FT,
                                 T::Array,
                               par::Array,
                             p_atm::FT,
                              p_O₂::FT,
                               r25::FT,
                         curvature::FT,
                                qy::FT) where {FT}
    # create lists to return
    N       = length(gsc)
    list_an = zeros(FT,N)
    list_ag = zeros(FT,N)
    list_r  = zeros(FT,N)
    list_pi = zeros(FT,N)

    # iterate through the list
    for i in 1:N
        _gsc = gsc[i]
        _par = par[i]
        _T   = T[i]

        # when there is no light
        if _par < 1e-3
            if r25==Inf
                r25 = model.VR * v25
            end
            tar_r  = get_r(model.ReT, r25, _T)
            tar_an = -tar_r
            tar_ag = FT(0.0)
            tar_pi = p_a + tar_r*FT(1e-6) * p_atm / _gsc
        else
            int_y  = p_a / p_atm * _gsc
            tar_pi = FT(0.1)
            tar_an = FT(0.0)
            tar_ag = FT(0.0)
            tar_r  = FT(0.0)
            count  = 0
            while true
                count  += 1
                an,ag,r = an_ag_r_from_pi(model,
                                          tar_pi,
                                          v25,
                                          p25,
                                          _par,
                                          _T,
                                          r25,
                                          qy)
                tar_g   = (int_y - an*FT(1e-6)) / tar_pi * p_atm

                if abs(tar_g - _gsc) < 1e-6
                    tar_an = an
                    tar_ag = ag
                    tar_r  = r
                    break
                end

                if count > 50
                    tar_an = an
                    tar_ag = ag
                    tar_r  = r
                    printstyled("Total number of iterations exceeds 50 times for C4ParaSet with Array!\n", color=:red)
                    println("\tThe provided PAR is ", par)
                    println("\tThe provided gsc is ", gsc)
                    println("\tThe provided V25 is ", v25)
                    println("")
                    break
                end

                tar_qi  = tar_pi + FT(0.001)
                ao,ag,r = an_ag_r_from_pi(model,
                                          tar_qi,
                                          v25,
                                          p25,
                                          _par,
                                          _T,
                                          r25,
                                          qy)
                tar_h   = (int_y - ao*FT(1e-6)) / tar_qi * p_atm

                slope   = (tar_h - tar_g) * 1000
                tar_pi += (_gsc - tar_g) / slope
            end
        end

        list_an[i] = tar_an
        list_ag[i] = tar_ag
        list_r[i]  = tar_r
        list_pi[i] = tar_pi
    end

    # return the an, ag, r, and p_i
    return list_an, list_ag, list_r, list_pi
end
