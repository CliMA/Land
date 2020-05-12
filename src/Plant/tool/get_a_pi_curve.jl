# this function is meant to obtain the Anet-Ci curve
# this function may also be used to test the phtosynthesis module
function get_a_pi_curve(v25=80.0, j25=135.0, gamma=2.5, tem=298.15, par=1200.0, p_o2=21278.25, r25=false, unit="K")
    # generate a list of p_i from gamma to 200 Pa
    list_pi = [val for val in gamma:1.0:200.0]
    list_ag = list_pi * -Inf
    list_an = list_pi * -Inf

    # iterate through the list
    for i in 1:length(list_pi)
        p_i        = list_pi[i]
        a_g        = get_leaf_ag_from_pi(v25,j25,gamma,p_i,tem,par,p_o2,r25,unit)
        a_n        = get_leaf_an_from_pi(v25,j25,gamma,p_i,tem,par,p_o2,r25,unit)
        list_ag[i] = a_g
        list_an[i] = a_n
    end

    # return the arrays
    return list_pi,list_an,list_ag
end
