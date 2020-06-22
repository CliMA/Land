#=
function get_p_base_q_list_from_q_rs(tree::Tree, flow::FT) where {FT}
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

This function use RootSolvers to get q_layer. A warning message will display if total iterations exceed 50 times.
"""
function get_q_layer_from_p_base_rs(root_layer::RootLayer{FT}, p_base::FT) where {FT}
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
# Calculate Weibull function considering temperature
#
###############################################################################

# TODO if p_his < p-crit, use 1e-12
function weibull_k_ratio(b::FT, c::FT, p::FT, T::FT, p_his::FT) where {FT}
    p_25  = min( p_his, p / relative_surface_tension(T) )
    kr_25 = max( FT(1e-6), exp( -1 * (-p_25/b) ^ (c) ) / relative_viscosity(T) )
    return kr_25
end

function weibull_k_leaf_ratio(leaf::Leaf{FT}, T::FT) where {FT}
    p_25  = leaf.p_element[end] / relative_surface_tension(T)
    kr_25 = max( FT(1e-6), exp( -1 * (-p_25 / leaf.b) ^ (leaf.c) ) / relative_viscosity(T) )
    return kr_25
end

=#
