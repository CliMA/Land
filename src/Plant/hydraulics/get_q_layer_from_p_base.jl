"""
    get_q_layer_from_p_base(root_layer, p_base)

Flow rate in a root layer `q_layer`, given
- `root_layer` One [`RootLayer`](@ref) in [`Root`](@ref) in [`Tree`](@ref)
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
