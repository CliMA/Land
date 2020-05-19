"""
    get_struct_from_q(root_layer, flow; p_ini)

End pressure `p_end` (in MPa) for [`RootLayer`](@ref), including impacts from rhizosphere conductance and gravity, given
- `root_layer` One [`RootLayer`](@ref) in [`Root`](@ref) in [`Tree`](@ref)
- `flow` Flow rate (in `mol s⁻¹`) in the given [`RootLayer`](@ref)
- `p_ini` Upstream soil water potential if `p_ini` is given, otherwise (`Inf`) `root_layer.p_ups` will be used
"""
function get_struct_p_end_from_q(root_layer::RootLayer, flow::FT; p_ini::FT=FT(Inf)) where {FT}
    if p_ini==Inf
        p_end = root_layer.p_ups
    else
        p_end = p_ini
    end

    #=
    compute the pressure drop along the rhizosphere, note that p_end is corrected for minor change in surface tension
    1. Using van Genuchten equation to compute soil hydraulic conductance Krhiz in the rhizosphere
    2. Calculate dP = flow / Krhiz_i for each of the 10 shells of the rhizosphere
    3. Update the P at the end of the rhizosphere
    =#
    a   = root_layer.soil_a
    m   = root_layer.soil_m
    n   = root_layer.soil_n
    k   = root_layer.k_rhiz
    k_s = get_relative_surface_tension(root_layer.t_soil)
    k_t = get_relative_viscosity(root_layer.t_soil)
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
    for i in 1:length(root_layer.p_element)
        p_25   = min(root_layer.p_history[i], p_end / get_relative_surface_tension(root_layer.t_element[i]))
        k_25   = root_layer.k_element[i] * exp( -1 * (-p_25/root_layer.b) ^ (root_layer.c) )
        k      = k_25 / get_relative_viscosity(root_layer.t_element[i])
        p_end -= flow / k + FT(ρ_H₂O) * FT(gravity) * root_layer.z_element[i] * FT(1e-6)
    end

    return p_end
end

"""
    get_struct_from_q(stem, flow; p_ini)

End pressure `p_end` (in MPa) for stem `Stem`, including impact from gravity, given
- `stem` [`Stem`](@ref) as Trunk in [`Tree`](@ref) or [`Stem`](@ref) as side branch in [`Branch`](@ref) in [`Tree`](@ref)
- `flow` Flow rate (in `mol s⁻¹`) in the given [`Stem`](@ref)
- `p_ini` Upstream xylem pressure if `p_ini` is given, otherwise (`Inf`) `stem.p_ups` will be used
"""
function get_struct_p_end_from_q(stem::Stem, flow::FT; p_ini::FT=(Inf)) where {FT}
    if p_ini==Inf
        p_end = stem.p_ups
    else
        p_end = p_ini
    end

    #=
    1. Compare the P with drought history in the xylem
    2. Using the Weibull function to compute the pressure drop along the xylem: k = k_max * exp( -(-P/b)^c )
    3. Calculate the pressure drop along each slice in the xylem, including the impact from gravity
    4. Update the P at the end of the xylem
    =#
    for i in 1:length(stem.p_element)
        p_25   = min(stem.p_history[i], p_end / get_relative_surface_tension(stem.t_element[i]))
        k_25   = stem.k_element[i] * exp( -1 * (-p_25/stem.b) ^ (stem.c) )
        k      = k_25 / get_relative_viscosity(stem.t_element[i])
        p_end -= flow / k + FT(ρ_H₂O) * FT(gravity) * stem.z_element[i] * FT(1e-6)
    end

    return p_end
end




"""
    get_struct_from_q(leaf, flow; p_ini)

End pressure `p_end` (in MPa) for leaf `leaf`, including impact from gravity, given
- `leaf` [`Leaf`](@ref) in [`CanopyLayer`](@ref) in [`Canopy`](@ref) in [`Tree`](@ref)
- `flow` Flow rate per leaf area (`mol m⁻² s⁻¹`) in the given [`Leaf`](@ref)
- `p_ini` Upstream xylem pressure if `p_ini` is given, otherwise (`Inf`) `leaf.p_ups` will be used
"""
function get_struct_p_end_from_q(leaf::Leaf, flow::FT; p_ini::FT=FT(Inf)) where {FT}
    if p_ini==Inf
        p_end = leaf.p_ups
    else
        p_end = p_ini
    end

    #=
    1. Compare the P with drought history in the xylem
    2. Using the Weibull function to compute the pressure drop along the xylem: k = k_max * exp( -(-P/b)^c )
    3. Calculate the pressure drop along each slice in the xylem, including the impact from gravity
    4. Update the P at the end of the xylem
    =#
    for i in 1:length(leaf.p_element)
        p_25   = min(leaf.p_history[i], p_end / get_relative_surface_tension(leaf.t_element[i]))
        k_25   = leaf.k_element[i] * exp( -1 * (-p_25/leaf.b) ^ (leaf.c) )
        k      = k_25 / get_relative_viscosity(leaf.t_element[i])
        p_end -= flow / k + FT(ρ_H₂O) * FT(gravity) * leaf.z_element[i] * FT(1e-6)
    end

    return p_end
end
