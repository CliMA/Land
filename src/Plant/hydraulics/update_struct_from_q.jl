"""
    update_struct_from_q!(root_layer, flow; p_ini)

Update root hydraulic information considering the impact from rhizosphere conductance and gravity, given
- `root_layer` One [`RootLayer`](@ref) in [`Root`](@ref) in [`Tree`](@ref)
- `flow` Flow rate in the given root layer
- `p_ini` Upstream soil water potential if `p_ini` is given, otherwise (`Inf`) `root_layer.p_ups` will be used
"""
function update_struct_from_q!(root_layer::RootLayer, flow::FT; p_ini::FT=FT(Inf)) where {FT}
    if p_ini==Inf
        p_end = root_layer.p_ups
    else
        p_end = p_ini
    end

    # update the flow rate in the struct
    root_layer.q = flow

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

    # update the water potential at rhizosphere-root connection
    root_layer.p_rhiz = p_end

    #=
    1. Compare the P with drought history in the xylem
    2. Using the Weibull function to compute the pressure drop along the xylem: k = k_max * exp( -(-P/b)^c )
    3. Calculate the pressure drop along each slice in the xylem, including the impact from gravity
    4. Update the P at the end of the xylem
    5. Also, update the drought history if the P is more negative than the previous history
    =#
    for i in 1:length(root_layer.p_element)
        p_25   = min(root_layer.p_history[i], p_end / get_relative_surface_tension(root_layer.t_element[i]))
        k_25   = root_layer.k_element[i] * exp( -1 * (-p_25/root_layer.b) ^ (root_layer.c) )
        k      = k_25 / get_relative_viscosity(root_layer.t_element[i])
        p_end -= flow / k + FT(ρ_H₂O) * FT(gravity) * root_layer.z_element[i] * FT(1e-6)

        # update the values in each element
        if p_25<root_layer.p_history[i]
            root_layer.p_history[i] = p_25
        end
        root_layer.p_element[i] = p_end
    end

    # update the xylem pressure at the tree base
    root_layer.p_dos = p_end
end




"""
    update_struct_from_q!(stem, flow; p_ini)

Update the stem hydraulic information considering the impact from gravity, given
- `stem` [`Stem`](@ref) as trunk in a [`Tree`](@ref) or as side branch in [`Branch`](ref) in [`Tree`](@ref)
- `flow::FT` Flow rate in the given stem
- `p_ini::FT` Upstream xylem pressure if `p_ini` is given, otherwise (`Inf`) `stem.p_ups` will be used
"""
function update_struct_from_q!(stem::Stem, flow::FT, p_ini::FT=FT(Inf)) where {FT}
    if p_ini==Inf
        p_end = stem.p_ups
    else
        p_end = p_ini
    end

    # update the flow rate in the struct
    stem.q = flow

    #=
    1. Compare the P with drought history in the xylem
    2. Using the Weibull function to compute the pressure drop along the xylem: k = k_max * exp( -(-P/b)^c )
    3. Calculate the pressure drop along each slice in the xylem, including the impact from gravity
    4. Update the P at the end of the xylem
    5. Also, update the drought history if the P is more negative than the previous history
    =#
    for i in 1:length(stem.p_element)
        p_25   = min(stem.p_history[i], p_end / get_relative_surface_tension(stem.t_element[i]))
        k_25   = stem.k_element[i] * exp( -1 * (-p_25/stem.b) ^ (stem.c) )
        k      = k_25 / get_relative_viscosity(stem.t_element[i])
        p_end -= flow / k + FT(ρ_H₂O) * FT(gravity) * stem.z_element[i] * FT(1e-6)

        # update the values in each element
        if p_25<stem.p_history[i]
            stem.p_history[i] = p_25
        end
        stem.p_element[i] = p_end
    end

    # update the xylem pressure at the tree base
    stem.p_dos = p_end
end




"""
    update_struct_from_q!(leaf, flow; p_ini)

Update the leaf hydraulic information considering the impact from gravity, given
- `leaf` [`Leaf`](@ref) type in [`CanopyLayer`](@ref) in a [`Tree`](@ref)
- `flow` Flow rate per leaf area in the given leaf
- `p_ini` Upstream xylem pressure if `p_ini` is given, otherwise (`Inf`) `leaf.p_ups` will be used
"""
function update_struct_from_q!(leaf::Leaf, flow::FT, p_ini::FT=FT(Inf)) where {FT}
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
    5. Also, update the drought history if the P is more negative than the previous history
    =#
    for i in 1:length(leaf.p_element)
        p_25   = min(leaf.p_history[i], p_end / get_relative_surface_tension(leaf.t_element[i]))
        k_25   = leaf.k_element[i] * exp( -1 * (-p_25/leaf.b) ^ (leaf.c) )
        k      = k_25 / get_relative_viscosity(leaf.t_element[i])
        p_end -= flow / k + FT(ρ_H₂O) * FT(gravity) * leaf.z_element[i] * FT(1e-6)

        # update the values in each element
        if p_25<leaf.p_history[i]
            leaf.p_history[i] = p_25
        end
        leaf.p_element[i] = p_end
    end

    # update the xylem pressure at the tree base
    leaf.p_dos = p_end
end
