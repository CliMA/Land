# this function calculates the q_layer from p_soil in the layer and p_base
# this function use newton raphson to get q_layer
function get_q_layer_from_p_base(root_layer::StructTreeRootLayer, p_base::Number)
    q_min   = 0.0
    q_max   = 1.0
    q_layer = 0.5

    # use while loop to get q_layer
    count   = 0
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
            q_layer  = 0.5 * (q_min + q_max)
        else
            p_1      = get_struct_p_end_from_q(root_layer, q_layer+1e-3)
            slope    = (p_1 - p_0) * 1e3
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