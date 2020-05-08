# this function returns the end pressure for a root, stem, or leaf struct
# flow in mol s^-1, p_end in MPa
function get_struct_p_end_from_q(any_struct, flow, p_ini=Inf)
    if p_ini>0
        p_end = any_struct.p_ini
    else
        p_end = p_ini
    end

    # compute the pressure drop along the xylem
    for i in 1:length(any_struct.p_element)
        p_25   = min(any_struct.p_history[i], p_end / get_relative_surface_tension(any_struct.t_element[i]))
        k_25   = any_struct.k_element[i] * exp( -1.0 * (-p_25/any_struct.b) ^ (any_struct.c) )
        k      = k_25 / get_relative_viscosity(any_struct.t_element[i])
        p_end -= flow / k + 998.0 * 9.8 * any_struct.h_element[i] * 1E-6
    end

    # return the result
    return p_end
end
