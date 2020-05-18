"""
    get_p_base_q_list_from_q(tree, flow)

# Arguments
- `tree::Tree`    Tree struct
- `flow::FT`      Flow rate for the whole tree, mol s⁻¹

# Description
This function returns the p_base from total transpiration rate using Newton Raphson.
A warning will display if total iterations is beyond 50 times.
"""
function get_p_base_q_list_from_q(tree::Tree, flow::FT) where {FT}
    p_base   = FT(0.0)
    q_list_0 = FT(0.0)

    # use while loop to get p_base, Δflow cannot be too small to avoid numerical issues
    Δflow = FT(1e-3)
    count = 0
    while true
        # calculate the q_list from an assumed p_base
        q_list_0 = [get_q_layer_from_p_base(root_layer, p_base) for root_layer in tree.roots.root_list]
        q_0      = sum(q_list_0)

        # break if sum of q equals flow rate
        if abs(q_0-flow)<1E-6
            break
        end

        # calculate the new q based on slope
        q_list_1 = [get_q_layer_from_p_base(root_layer, p_base+Δflow) for root_layer in tree.roots.root_list]
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
