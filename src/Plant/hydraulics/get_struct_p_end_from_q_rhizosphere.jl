# this function is specifically for root with rhizosphere conductance
# flow in mol s^-1, p_end in MPa
function get_struct_p_end_from_q_rhizosphere(any_struct::struct_tree_root, flow, p_ini=Inf)
    if p_ini>0
        p_end = any_struct.p_ups
    else
        p_end = p_ini
    end

    # compute the pressure drop along the rhizosphere, note that p_end is corrected for minor change in surface tension
    a   = any_struct.soil_a
    m   = any_struct.soil_m
    n   = any_struct.soil_n
    k   = any_struct.k_rhiz
    k_s = get_relative_surface_tension(any_struct.t_soil)
    k_t = get_relative_viscosity(any_struct.t_soil)
    for i in 0:9
        if p_end<0
            shell_t = (1.0 / (1.0 + (a*(-p_end/k_s))^n)) ^ m
        else
            shell_t = 1.0
        end
        shell_f  = sqrt(shell_t) * (1.0 - (1.0-shell_t^(1.0/m)) ^ m) ^ 2.0
        shell_k  = k * shell_f * log(10.0) / log((10.0-0.9*i)/(10.0-0.9*(i+1))) / k_t
        dp       = flow / shell_k
        p_end   -= dp
    end

    # compute the pressure drop along the xylem
    for i in 1:length(any_struct.p_element)
        p_25   = min(any_struct.p_history[i], p_end / get_relative_surface_tension(any_struct.t_element[i]))
        k_25   = any_struct.k_element[i] * exp( -1.0 * (-p_25/any_struct.b) ^ (any_struct.c) )
        k      = k_25 / get_relative_viscosity(any_struct.t_element[i])
        p_end -= flow / k + 998.0 * 9.8 * any_struct.z_element[i] * 1E-6
    end

    # return the result
    return p_end
end
