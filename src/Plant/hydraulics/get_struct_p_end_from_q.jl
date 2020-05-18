"""
    get_struct_from_q(any_struct, flow; p_ini)

# Arguments
- `any_struct::RootLayer`    One root layer in Tree struct
- `flow::FT`                 Flow rate in the given root layer
- `p_ini::FT`                Upstream soil water potential if p_ini is given, otherwise (Inf) any_struct.p_ups will be used

# Description
This function returns the end pressure for rootlayer struct, including impacts from rhizosphere conductance and gravity.
Flow in mol s⁻¹, p_end in MPa
"""
function get_struct_p_end_from_q(any_struct::RootLayer, flow::FT; p_ini::FT=FT(Inf)) where {FT}
    if p_ini==Inf
        p_end = any_struct.p_ups
    else
        p_end = p_ini
    end

    #=
    compute the pressure drop along the rhizosphere, note that p_end is corrected for minor change in surface tension
    1. Using van Genuchten equation to compute soil hydraulic conductance Krhiz in the rhizosphere
    2. Calculate dP = flow / Krhiz_i for each of the 10 shells of the rhizosphere
    3. Update the P at the end of the rhizosphere
    =#
    a   = any_struct.soil_a
    m   = any_struct.soil_m
    n   = any_struct.soil_n
    k   = any_struct.k_rhiz
    k_s = get_relative_surface_tension(any_struct.t_soil)
    k_t = get_relative_viscosity(any_struct.t_soil)
    for i in 0:9
        if p_end<=0
            shell_t = (1 / (1 + (a*(-p_end/k_s))^n)) ^ m
        else
            shell_t = FT(1.0)
        end
        shell_f  = sqrt(shell_t) * (1 - (1-shell_t^(1/m)) ^ m) ^ 2
        shell_k  = k * shell_f * log( FT(10) ) / log( FT((10-0.9*i) / (10-0.9*(i+1))) ) / k_t
        dp       = flow / shell_k
        p_end   -= dp
    end

    #=
    1. Compare the P with drought history in the xylem
    2. Using the Weibull function to compute the pressure drop along the xylem: k = k_max * exp( -(-P/b)^c )
    3. Calculate the pressure drop along each slice in the xylem, including the impact from gravity
    4. Update the P at the end of the xylem
    =#
    for i in 1:length(any_struct.p_element)
        p_25   = min(any_struct.p_history[i], p_end / get_relative_surface_tension(any_struct.t_element[i]))
        k_25   = any_struct.k_element[i] * exp( -1 * (-p_25/any_struct.b) ^ (any_struct.c) )
        k      = k_25 / get_relative_viscosity(any_struct.t_element[i])
        p_end -= flow / k + ρ_H₂O * gravity * any_struct.z_element[i] * FT(1e-6)
    end

    # return the result
    return p_end
end




"""
    get_struct_from_q(any_struct, flow; p_ini)

# Arguments
- `any_struct::Stem`    Trunk or Branch in Tree struct
- `flow::FT`            Flow rate in the given Stem
- `p_ini::FT`           Upstream xylem pressure if p_ini is given, otherwise (Inf) any_struct.p_ups will be used
    
# Description
This function returns the end pressure for Stem struct, including impact from gravity.
Flow in mol s⁻¹, p_end in MPa
"""
function get_struct_p_end_from_q(any_struct::Stem, flow::FT; p_ini::FT=(Inf)) where {FT}
    if p_ini==Inf
        p_end = any_struct.p_ups
    else
        p_end = p_ini
    end

    #=
    1. Compare the P with drought history in the xylem
    2. Using the Weibull function to compute the pressure drop along the xylem: k = k_max * exp( -(-P/b)^c )
    3. Calculate the pressure drop along each slice in the xylem, including the impact from gravity
    4. Update the P at the end of the xylem
    =#
    for i in 1:length(any_struct.p_element)
        p_25   = min(any_struct.p_history[i], p_end / get_relative_surface_tension(any_struct.t_element[i]))
        k_25   = any_struct.k_element[i] * exp( -1 * (-p_25/any_struct.b) ^ (any_struct.c) )
        k      = k_25 / get_relative_viscosity(any_struct.t_element[i])
        p_end -= flow / k + ρ_H₂O * gravity * any_struct.z_element[i] * FT(1e-6)
    end

    # return the result
    return p_end
end




"""
    get_struct_from_q(any_struct, flow; p_ini)

# Arguments
- `any_struct::Leaf`    Leaf in CanopyLayer in a Tree struct
- `flow::FT`            Flow rate in the given Leaf
- `p_ini::FT`           Upstream xylem pressure if p_ini is given, otherwise (Inf) any_struct.p_ups will be used
    
# Description
This function returns the end pressure for Leaf struct, including impact from gravity.
Flow in mol s⁻¹, p_end in MPa
"""
function get_struct_p_end_from_q(any_struct::Leaf, flow::FT; p_ini::FT=FT(Inf)) where {FT}
    if p_ini==Inf
        p_end = any_struct.p_ups
    else
        p_end = p_ini
    end

    #=
    1. Compare the P with drought history in the xylem
    2. Using the Weibull function to compute the pressure drop along the xylem: k = k_max * exp( -(-P/b)^c )
    3. Calculate the pressure drop along each slice in the xylem, including the impact from gravity
    4. Update the P at the end of the xylem
    =#
    for i in 1:length(any_struct.p_element)
        p_25   = min(any_struct.p_history[i], p_end / get_relative_surface_tension(any_struct.t_element[i]))
        k_25   = any_struct.k_element[i] * exp( -1 * (-p_25/any_struct.b) ^ (any_struct.c) )
        k      = k_25 / get_relative_viscosity(any_struct.t_element[i])
        p_end -= flow / k + ρ_H₂O * gravity * any_struct.z_element[i] * FT(1e-6)
    end

    # return the result
    return p_end
end
