# function to get marginal carbon penalty for the Wang 2020 model
function get_marginal_penalty_wang(canopyi::StructTreeCanopyLayer, indx::Number)
    ∂Θ∂E = canopyi.an_list[indx] / (canopyi.ec_list[indx] - canopyi.e_list[indx]) * 1e-6
    return ∂Θ∂E
end




# function to get marginal carbon penalty list for the Wang 2020 model
function get_marginal_penalty_wang(canopyi::StructTreeCanopyLayer)
    list_∂Θ∂E = canopyi.an_list ./ (canopyi.ec_list .- canopyi.e_list) .* 1e-6
    list_∂Θ∂E[ list_∂Θ∂E .<=0 ] .= Inf
    return list_∂Θ∂E
end