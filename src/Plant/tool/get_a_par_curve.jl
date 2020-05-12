# this function is meant to obtain the Anet-PAR curve
# this function may also be used to test the phtosynthesis module
function get_a_par_curve(v25=80.0, j25=135.0, gamma=2.5, gsc=0.2, ca=40.0, tem=298.15, p_atm=101325.0,p_o2=21278.25, r25=false, unit="K")
    # generate a list of p_i from gamma to 200 Pa
    list_par = [val for val in 5.0:5.0:2000.0]
    list_ag  = list_par * -Inf
    list_an  = list_par * -Inf
    list_pi  = list_par * -Inf

    # iterate through the list
    for i in 1:length(list_par)
        par         = list_par[i]
        a_n,p_i,a_g = get_leaf_an_pi_ag_from_gsc(v25,j25,gamma,gsc,ca,tem,par,p_atm,p_o2,r25,unit)
        list_ag[i]  = a_g
        list_an[i]  = a_n
        list_pi[i]  = p_i
    end

    # return the arrays
    return list_par,list_an,list_ag,list_pi
end
