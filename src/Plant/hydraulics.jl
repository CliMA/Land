###############################################################################
#
# Calculate P_base and list for each root layer from a known Q
# These functions passed the FT test
# These functions are documented in the Plant page
#
###############################################################################
"""
    get_p_base_q_list_from_q(tree, flow)

Tree base pressyre `p_base` and a list of flow rate in each root layer `q_list`, given
- `tree` A [`Tree`](@ref) type
- `flow` Flow rate for the whole tree `[mol s⁻¹]`

This function returns the p_base from total transpiration rate using Newton Raphson. A warning will display if total iterations is beyond 50 times.
"""
function get_p_base_q_list_from_q(tree::Tree, flow::FT) where {FT}
    p_base   = FT(0.0)
    q_list_0 = FT(0.0)

    # use while loop to get p_base, Δflow cannot be too small to avoid numerical issues
    Δflow = FT(1e-3)
    count = 0
    while true
        # calculate the q_list from an assumed p_base
        q_list_0 = [get_q_layer_from_p_base(root_layer, p_base) for root_layer in tree.root_list]
        q_0      = sum(q_list_0)

        # break if sum of q equals flow rate
        if abs(q_0-flow)<1E-6
            break
        end

        # calculate the new q based on slope
        q_list_1 = [get_q_layer_from_p_base(root_layer, p_base+Δflow) for root_layer in tree.root_list]
        q_1      = sum(q_list_1)
        slope    = (q_1-q_0) / Δflow
        p_base  += (flow-q_0) / slope

        # break if total iterations >= 50
        count += 1
        if count>=50
            println("Warning: total number of iteration exceeds 50 times when computing p_base and flow rates in the root layers from a total flow rate, please check whether the result is reliable!")
            println("The model predicted total flow rate is ", q_0)
            println("The given total flow rate is ", flow)
            break
        end
    end

    return p_base,q_list_0
end




"""
    get_q_layer_from_p_base(root_layer, p_base)

Flow rate in a root layer `q_layer`, given
- `root_layer` One [`RootLayer`](@ref) in [`Tree`](@ref)
- `p_base` Xylem water pressure at the tree base

This function use Newton Raphson combined with Bi-section to get q_layer. A warning message will display if total iterations exceed 50 times.
"""
function get_q_layer_from_p_base(root_layer::RootLayer, p_base::FT) where {FT}
    q_min   = FT(0.0)
    q_max   = FT(1.0)
    q_layer = FT(0.5)

    # use while loop to get q_layer, Δp cannot be too small to avoid numerical issues
    Δp    = FT(1e-3)
    count = 0
    while true
        p_0 = get_struct_p_end_from_q(root_layer, q_layer)

        # if meet the requirement break
        if abs(p_0 - p_base) < 1E-6
            break
        end

        # update q_min and q_max
        if p_0 < p_base
            q_max = q_layer
        else
            q_min = q_layer
        end
        
        # update q_layer from bi-section or newton raphson
        if p_0 <= -20.0
            q_layer  = FT(0.5) * (q_min + q_max)
        else
            p_1      = get_struct_p_end_from_q(root_layer, q_layer+Δp)
            slope    = (p_1 - p_0) / Δp
            q_layer += (p_base - p_0) / slope
        end

        # break if total iterations >= 50 or branch into bi-section method?
        count += 1
        if count>=50
            println("Warning: total number of iteration exceeds 50 times when computing flow rate from p_base, please check whether the result is reliable!")
            println("The model predicted p_base is ", p_0)
            println("The given p_base is ", p_base)
            break
        end
    end

    return q_layer
end








###############################################################################
#
# Calculate or Update Pressure profile along a tree hydraulic element
# These functions passed the FT test
# These functions are documented in the Plant page
#
###############################################################################
"""
    get_struct_p_end_from_q(root_layer::RootLayer, flow::FT; p_ini::FT=FT(Inf))
    get_struct_p_end_from_q(stem::Stem, flow::FT; p_ini::FT=FT(Inf))
    get_struct_p_end_from_q(leaf::Leaf, flow::FT; p_ini::FT=FT(Inf))

End pressure `p_end` (in MPa), including impacts from rhizosphere conductance and gravity, given
- `root_layer` One [`RootLayer`](@ref) in [`Tree`](@ref)
- `stem` [`Stem`](@ref) as Trunk in [`Tree`](@ref) or [`Stem`](@ref) as side branch in [`Tree`](@ref)
- `leaf` [`Leaf`](@ref) in [`CanopyLayer`](@ref) in [`Tree`](@ref)
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
    k_s = relative_surface_tension(root_layer.t_soil)
    k_t = relative_viscosity(root_layer.t_soil)
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
        p_25   = min(root_layer.p_history[i], p_end / relative_surface_tension(root_layer.t_element[i]))
        k_25   = root_layer.k_element[i] * exp( -1 * (-p_25/root_layer.b) ^ (root_layer.c) )
        k      = k_25 / relative_viscosity(root_layer.t_element[i])
        p_end -= flow / k + FT(ρ_H₂O) * FT(GRAVITY) * root_layer.z_element[i] * FT(1e-6)
    end

    return p_end
end

function get_struct_p_end_from_q(stem::Stem, flow::FT; p_ini::FT=FT(Inf)) where {FT}
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
        p_25   = min(stem.p_history[i], p_end / relative_surface_tension(stem.t_element[i]))
        k_25   = stem.k_element[i] * exp( -1 * (-p_25/stem.b) ^ (stem.c) )
        k      = k_25 / relative_viscosity(stem.t_element[i])
        p_end -= flow / k + FT(ρ_H₂O) * FT(GRAVITY) * stem.z_element[i] * FT(1e-6)
    end

    return p_end
end

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
        p_25   = min(leaf.p_history[i], p_end / relative_surface_tension(leaf.t_element[i]))
        k_25   = leaf.k_element[i] * exp( -1 * (-p_25/leaf.b) ^ (leaf.c) )
        k      = k_25 / relative_viscosity(leaf.t_element[i])
        p_end -= flow / k + FT(ρ_H₂O) * FT(GRAVITY) * leaf.z_element[i] * FT(1e-6)
    end

    return p_end
end




"""
    update_struct_from_q!(root_layer::RootLayer, flow::FT; p_ini::FT=FT(Inf))
    update_struct_from_q!(stem::Stem, flow::FT, p_ini::FT=FT(Inf))
    update_struct_from_q!(leaf::Leaf, flow::FT, p_ini::FT=FT(Inf))

Update hydraulic information considering the impact from rhizosphere conductance and gravity, given
- `root_layer` One [`RootLayer`](@ref) in [`Tree`](@ref)
- `stem` [`Stem`](@ref) as trunk in a [`Tree`](@ref) or as side branch in [`Tree`](@ref)
- `leaf` [`Leaf`](@ref) type in [`CanopyLayer`](@ref) in a [`Tree`](@ref)
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
    k_s = relative_surface_tension(root_layer.t_soil)
    k_t = relative_viscosity(root_layer.t_soil)
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
        p_25   = min(root_layer.p_history[i], p_end / relative_surface_tension(root_layer.t_element[i]))
        k_25   = root_layer.k_element[i] * exp( -1 * (-p_25/root_layer.b) ^ (root_layer.c) )
        k      = k_25 / relative_viscosity(root_layer.t_element[i])
        p_end -= flow / k + FT(ρ_H₂O) * FT(GRAVITY) * root_layer.z_element[i] * FT(1e-6)

        # update the values in each element
        if p_25<root_layer.p_history[i]
            root_layer.p_history[i] = p_25
        end
        root_layer.p_element[i] = p_end
    end

    # update the xylem pressure at the tree base
    root_layer.p_dos = p_end
end

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
        p_25   = min(stem.p_history[i], p_end / relative_surface_tension(stem.t_element[i]))
        k_25   = stem.k_element[i] * exp( -1 * (-p_25/stem.b) ^ (stem.c) )
        k      = k_25 / relative_viscosity(stem.t_element[i])
        p_end -= flow / k + FT(ρ_H₂O) * FT(GRAVITY) * stem.z_element[i] * FT(1e-6)

        # update the values in each element
        if p_25<stem.p_history[i]
            stem.p_history[i] = p_25
        end
        stem.p_element[i] = p_end
    end

    # update the xylem pressure at the tree base
    stem.p_dos = p_end
end

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
        p_25   = min(leaf.p_history[i], p_end / relative_surface_tension(leaf.t_element[i]))
        k_25   = leaf.k_element[i] * exp( -1 * (-p_25/leaf.b) ^ (leaf.c) )
        k      = k_25 / relative_viscosity(leaf.t_element[i])
        p_end -= flow / k + FT(ρ_H₂O) * FT(GRAVITY) * leaf.z_element[i] * FT(1e-6)

        # update the values in each element
        if p_25<leaf.p_history[i]
            leaf.p_history[i] = p_25
        end
        leaf.p_element[i] = p_end
    end

    # update the xylem pressure at the tree base
    leaf.p_dos = p_end
end







###############################################################################
#
# Update tree e_crit information for each leaf
# This function passed the FT test
# This function is documented in the Plant page
#
###############################################################################
"""
    update_tree_e_crit!(tree::Tree; displaying::Bool=false)

Update the critical transpiration rate information, given
- `tree` A [`Tree`](@ref) type
- `displaying` If true, messages of how `e_crit` evolves step-by-step will display

This function update `e_crit` for each layer (sunlit and shaded). The `e_crit` for each leaf is either increased or decreased so that `p_leaf` is in the range -Inf to -20.0 MPa. The `e_crit increases by 1e-5 mol m⁻² s⁻¹` if `p_leaf` is less negative than -20 MPa. The `e_crit decreases by 1e-6 mol m⁻² s⁻¹` if `p_leaf` is -Inf.
"""
function update_tree_e_crit!(tree::Tree; displaying::Bool=false)
    # unpack necessary structs
    @unpack branch_list,canopy_list,root_list   = tree

    while true
        # define judge = 0
        judge = 0

        # q for the whole tree is the sum of q in each leaves
        q_branch_list = [sum(canopyi.ec_list .* canopyi.la_list) for canopyi in canopy_list]
        q_sum = sum( q_branch_list )

        # calculate the p_base from q_sum
        p_base,q_list = get_p_base_q_list_from_q(tree, q_sum)

        # calculate the trunk-branch joint pressure from q_sum and p_base
        p_trunk_branch = get_struct_p_end_from_q(tree.trunk, q_sum; p_ini=p_base)

        # for each canopy layer
        for indx in 1:length(branch_list)
            branchi = branch_list[indx]
            canopyi = canopy_list[indx]
            # compute branch leaf joint pressure for each canopy layer
            p_branch_leaf = get_struct_p_end_from_q(branchi, q_branch_list[indx]; p_ini=p_trunk_branch)

            # compute p_leaf for each leaf
            for indy in 1:length(canopyi.leaf_list)
                leaf = canopyi.leaf_list[indy]
                p_leaf = get_struct_p_end_from_q(leaf, canopyi.ec_list[indy]; p_ini=p_branch_leaf)
                # increas e_crit by 1E-6 or decrease it by 1E-7
                if p_leaf >= -20.0
                    judge += 1
                    canopyi.ec_list[indy] += 1e-5
                elseif p_leaf == -Inf
                    judge += 1
                    canopyi.ec_list[indy] -= 1e-6
                end
            end
        end

        # if judge > 0, continue
        if judge==0
            break
        else
            if displaying
                println("The Q_sum was ", q_sum, " and then ", judge, " modification(s) was(were) made.")
            end
        end
    end
end








###############################################################################
#
# Update tree e_crit information for each leaf
# This function passed the FT test
# This function is documented in the Plant page
#
###############################################################################
"""
    update_leaf_ak_max!(tree::Tree)

Update leaf maximal `A_net` allowed at the given environmental conditions, given
-`tree` A [`Tree`](@ref) type struct

This function is meant to work with Sperry et al. (2017) model, which requires a_max in the ∂Θ∂E function.

"""
function update_leaf_ak_max!(tree::Tree)
    @unpack canopy_list,photo_para_set = tree
    FTYP = eltype(tree.ba)

    # 1. iterate through canopy layers
    for canopyi in canopy_list
        # 2. iterate through the leaves in the canopy layers
        for i in 1:length(canopyi.leaf_list)
            leaf = canopyi.leaf_list[i]
            
            # compute a_max
            e_c  = canopyi.ec_list[i]
            gsw  = min( e_c / (canopyi.d_list[i] / canopyi.p_atm), canopyi.g_max )
            gsc  = gsw / FTYP(1.6) / ( 1 + canopyi.g_ias_c * gsw^(canopyi.g_ias_e) )
            re_c = an_ag_r_pi_from_gsc(
                                       photo_para_set,
                                       gsc,
                                       canopyi.v_max,
                                       canopyi.j_max,
                                       canopyi.p_max,
                                       canopyi.p_a,
                                       canopyi.t_list[i],
                                       canopyi.par_list[i],
                                       canopyi.p_atm,
                                       canopyi.p_O₂,
                                       canopyi.r_25,
                                       canopyi.curvature,
                                       canopyi.qy)
            a_m  = re_c[1]

            # compute k_max
            p_m  = get_struct_p_end_from_q(leaf, FTYP(0.0))
            k_m  = exp( -(-p_m/leaf.b)^leaf.c )

            # update a_max and k_max
            canopyi.am_list[i] = a_m
            canopyi.km_list[i] = k_m
        end
    end
end
