# function to get marginal carbon penalty for the Wang 2020 model
function get_marginal_penalty_wang(leaf::StructTreeLeaf)
    ∂Θ∂E = leaf.a_net / (leaf.e_crit - leaf.e) * 1e-6
    return ∂Θ∂E
end