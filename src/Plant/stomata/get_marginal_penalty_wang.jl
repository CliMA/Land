"""
    get_marginal_penalty_wang(canopyi, indx)

# Arguments
- `canopyi::CanopyLayer`    A canopy layer in the Tree Struct
- `indx::Number`            indx'th leaf in canopyi.leaf_list

# Description
This function returns marginal carbon penalty for the indx th leaf in the canopyi using the Wang 2020 model.
∂Θ∂E = A / (e_crit - e_leaf), unit mol mol⁻¹
Should be updated to Abstract Type once more models are added
"""
function get_marginal_penalty_wang(canopyi::CanopyLayer, indx::Number)
    FTYP = eltype(canopyi.an_list[indx])
    ∂Θ∂E = canopyi.an_list[indx] / (canopyi.ec_list[indx] - canopyi.e_list[indx]) * FTYP(1e-6)
    return ∂Θ∂E
end




"""
    get_marginal_penalty_wang(canopyi, indx)

# Arguments
- `canopyi::CanopyLayer`    A canopy layer in the Tree Struct

# Description
This function returns marginal carbon penalty for all the leavs in canopyi using the Wang 2020 model.
∂Θ∂E = A / (e_crit - e_leaf), unit mol mol⁻¹
Should be updated to Abstract Type once more models are added
"""
function get_marginal_penalty_wang(canopyi::CanopyLayer)
    FTYP      = eltype(canopyi.an_list[1])
    list_∂Θ∂E = canopyi.an_list ./ (canopyi.ec_list .- canopyi.e_list) .* FTYP(1e-6)
    list_∂Θ∂E[ list_∂Θ∂E .<= 0 ] .= FTYP(Inf)
    return list_∂Θ∂E
end