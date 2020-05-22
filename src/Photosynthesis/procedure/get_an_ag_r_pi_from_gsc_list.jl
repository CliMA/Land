"""
    get_an_ag_r_pi_from_gsc_list(; model, v25, j25, gsc_list, p_a, t_leaf_list, par_list, p_atm, p_O₂, r25, curvature, qy)

List of net photosynthetic rate `list_an`, gross photosynthetic rate `list_ag`, respiration rate `list_re`, and leaf internal CO₂ partial pressure `list_pi`, given
- `model` A C3ParaSet type parameter set
- `v25` Maixmal carboxylation rate at 298.15 K (25 Celcius)
- `j25` Maximal electron transport rate at 298.15 K (25 Celcius)
- `gsc_list` List of leaf diffusive conductance to CO₂
- `p_a` Atmospheric CO₂ partial pressure
- `t_leaf_list` List of leaf temperature
- `par_list` List of photosynthetic active radiation
- `p_atm` Atmospheric pressure
- `p_O₂` O₂ partial pressure
- `r25` Leaf respiration rate at 298.15 K (25 Celcius). If ==Inf, r25 will be computed from v25
- `curvature` Curvature parameter to calculate J
- `qy` Quantum yield of electron


"""
function get_an_ag_r_pi_from_gsc_list(;
                                      model::C3ParaSet   = C3VcJBernacchi{FT}(),
                                        v25::FT          = FT(80.0),
                                        j25::FT          = FT(135.0),
                                   gsc_list::Array{FT,1} = FT(0.1) .* ones(FT,10),
                                        p_a::FT          = FT(40.0),
                                t_leaf_list::Array{FT,1} = FT(298.15).* ones(FT,10),
                                   par_list::Array{FT,1} = FT(1000.0).* ones(FT,10),
                                      p_atm::FT          = FT(101325.0),
                                       p_O₂::FT          = FT(21278.25),
                                        r25::FT          = FT(Inf),
                                  curvature::FT          = FT(0.9),
                                         qy::FT          = FT(0.4081632653061224)) where {FT}
    # define lists of results
    len     = length(gsc_list)
    list_an = zeros(FT,len)
    list_ag = zeros(FT,len)
    list_re = zeros(FT,len)
    list_pi = zeros(FT,len)

    # iterate the list
    for indx in 1:length( len )
        gsc         =    gsc_list[indx]
        t_leaf      = t_leaf_list[indx]
        par         =    par_list[indx]
        an,ag,r,p_i = get_an_ag_r_pi_from_gsc(
                                              model = model,
                                                v25 = v25,
                                                j25 = j25,
                                                gsc = gsc,
                                                p_a = p_a,
                                             t_leaf = t_leaf,
                                                par = par,
                                              p_atm = p_atm,
                                               p_O₂ = p_O₂,
                                                r25 = r25,
                                          curvature = curvature,
                                                 qy = qy)
        list_an[indx] = an
        list_ag[indx] = ag
        list_re[indx] = r
        list_pi[indx] = p_i
    end

    # return the lists
    return list_an, list_ag, list_re, list_pi
end




"""
    get_an_ag_r_pi_from_gsc_list(; model, v25, p25, gsc, p_a, t_leaf, par, p_atm, r25, qy)

List of net photosynthetic rate `list_an`, gross photosynthetic rate `list_ag`, respiration rate `list_re`, and leaf internal CO₂ partial pressure `list_pi`, given
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

For C4 plants.

"""
function get_an_ag_r_pi_from_gsc_list(;
                                      model::C4ParaSet   = C4VcVpJCLM{FT}(),
                                        v25::FT          = FT(80.0),
                                        p25::FT          = FT(120.0),
                                   gsc_list::Array{FT,1} = FT(0.1).* ones(FT,10),
                                        p_a::FT          = FT(40.0),
                                t_leaf_list::Array{FT,1} = FT(298.15).* ones(FT,10),
                                   par_list::Array{FT,1} = FT(1000.0).* ones(FT,10),
                                      p_atm::FT          = FT(101325.0),
                                        r25::FT          = FT(Inf),
                                         qy::FT          = FT(0.4081632653061224)) where {FT}
    # define lists of results
    len     = length(gsc_list)
    list_an = zeros(FT,len)
    list_ag = zeros(FT,len)
    list_re = zeros(FT,len)
    list_pi = zeros(FT,len)

    # iterate the list
    for indx in 1:length( len )
        gsc         =    gsc_list[indx]
        t_leaf      = t_leaf_list[indx]
        par         =    par_list[indx]
        an,ag,r,p_i = get_an_ag_r_pi_from_gsc(
                                              model = model,
                                                v25 = v25,
                                                p25 = p25,
                                                gsc = gsc,
                                                p_a = p_a,
                                             t_leaf = t_leaf,
                                                par = par,
                                              p_atm = p_atm,
                                                r25 = r25,
                                                 qy = qy)
        list_an[indx] = an
        list_ag[indx] = ag
        list_re[indx] = r
        list_pi[indx] = p_i
    end

    # return the lists
    return list_an, list_ag, list_re, list_pi
end
