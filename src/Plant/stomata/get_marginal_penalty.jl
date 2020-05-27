"""
    get_marginal_penalty(canopyi, indx, scheme::OSMEller)

Marginal water penalty (`∂Θ/∂E`) for opening the stomata for a leaf in canopy layer, given
- `canopyi` A [`CanopyLayer`](@ref) in [`Canopy`](@ref) in a [`Tree`](@ref)
- `indx` indx'th [`Leaf`](@ref) in the `canopyi.leaf_list`
- `scheme` optimizaiton stomatal model scheme OSMEller

This function uses the Eller et al. (2018) model. `∂Θ/∂E = A / K * dK/dE`, unit `[mol mol⁻¹]`.
To make the Eller model suitable for multiple layer canopy, K at the leaf xylem end is used. Also, to make the model capable of handling drought history, loss of conductance is evaluated assuming no drought legacy from a more negative `p_leaf`.

"""
function get_marginal_penalty(canopyi::CanopyLayer, indx::Number, scheme::OSMEller)
    leaf = canopyi.leaf_list[indx]
    FTYP = eltype(leaf.b)

    # calculate the e, p, and k
    de   = FTYP(1e-6)
    e_0  = canopyi.e_list[indx]
    e_1  = e_0 + de
    p_0  = leaf.p_element[end]
    p_1  = get_struct_p_end_from_q(leaf, e_1)
    k_0  = exp( -(-p_0/leaf.b)^leaf.c )
    k_1  = exp( -(-p_1/leaf.b)^leaf.c )
    ∂Θ∂E = (k_0 - k_1) / de * canopyi.an_list[indx] / leaf.k_element[indx] * leaf.k_sla * 10
    return ∂Θ∂E
end




"""
    get_marginal_penalty(canopyi, scheme::OSMEller)

List of marginal water penalty for all the leaves, given
- `canopyi` A [`CanopyLayer`](@ref) in [`Canopy`](@ref) in a [`Tree`](@ref)
- `scheme` optimizaiton stomatal model scheme OSMEller

This function uses the Eller et al. (2018) model. `∂Θ/∂E = A / K * dK/dE`, unit `[mol mol⁻¹]`.
To make the Eller model suitable for multiple layer canopy, K at the leaf xylem end is used. Also, to make the model capable of handling drought history, loss of conductance is evaluated assuming no drought legacy from a more negative `p_leaf`.

"""
function get_marginal_penalty(canopyi::CanopyLayer, scheme::OSMEller)
    FTYP      = eltype(canopyi.la)
    list_∂Θ∂E = canopyi.an_list .* 0
    for indx in 1:length(canopyi.leaf_list)
        list_∂Θ∂E[indx] = get_marginal_penalty(canopyi, indx, scheme)
    end
    list_∂Θ∂E[ list_∂Θ∂E .<= 0 ] .= FTYP(Inf)
    return list_∂Θ∂E
end




"""
    get_marginal_penalty(canopyi, indx, scheme::OSMSperry)

Marginal water penalty (`∂Θ/∂E`) for opening the stomata for a leaf in canopy layer, given
- `canopyi` A [`CanopyLayer`](@ref) in [`Canopy`](@ref) in a [`Tree`](@ref)
- `indx` indx'th [`Leaf`](@ref) in the `canopyi.leaf_list`
- `scheme` optimizaiton stomatal model scheme OSMSperry

This function uses the Sperry et al. (2017) model. `∂Θ/∂E = A_max / K_max * dK/dE`, unit `[mol mol⁻¹]`.
To make the Sperry model suitable for multiple layer canopy, K at the leaf xylem end is used. Also, to make the model capable of handling drought history, loss of conductance is evaluated assuming no drought legacy from a more negative `p_leaf`.

"""
function get_marginal_penalty(canopyi::CanopyLayer, indx::Number, scheme::OSMSperry)
    leaf = canopyi.leaf_list[indx]
    FTYP = eltype(leaf.b)

    # calculate the e, p, and k
    de   = FTYP(1e-6)
    e_0  = canopyi.e_list[indx]
    e_1  = e_0 + de
    p_0  = leaf.p_element[end]
    p_1  = get_struct_p_end_from_q(leaf, e_1)
    k_0  = exp( -(-p_0/leaf.b)^leaf.c )
    k_1  = exp( -(-p_1/leaf.b)^leaf.c )
    ∂Θ∂E = (k_0 - k_1) / de * canopyi.am_list[indx] / canopyi.km_list[indx]
    return ∂Θ∂E
end




"""
    get_marginal_penalty(canopyi, scheme::OSMSperry)

List of marginal water penalty for all the leaves, given
- `canopyi` A [`CanopyLayer`](@ref) in [`Canopy`](@ref) in a [`Tree`](@ref)
- `scheme` optimizaiton stomatal model scheme OSMSperry

This function uses the Sperry et al. (2017) model. `∂Θ/∂E = A_max / K_max * dK/dE`, unit `[mol mol⁻¹]`.
To make the Sperry model suitable for multiple layer canopy, K at the leaf xylem end is used. Also, to make the model capable of handling drought history, loss of conductance is evaluated assuming no drought legacy from a more negative `p_leaf`.

"""
function get_marginal_penalty(canopyi::CanopyLayer, scheme::OSMSperry)
    FTYP      = eltype(canopyi.la)
    list_∂Θ∂E = canopyi.an_list .* 0
    for indx in 1:length(canopyi.leaf_list)
        list_∂Θ∂E[indx] = get_marginal_penalty(canopyi, indx, scheme)
    end
    list_∂Θ∂E[ list_∂Θ∂E .<= 0 ] .= FTYP(Inf)
    return list_∂Θ∂E
end




"""
    get_marginal_penalty(canopyi, indx, scheme::OSMWang)

Marginal water penalty (`∂Θ/∂E`) for opening the stomata for a leaf in canopy layer, given
- `canopyi` A [`CanopyLayer`](@ref) in [`Canopy`](@ref) in a [`Tree`](@ref)
- `indx` indx'th [`Leaf`](@ref) in the `canopyi.leaf_list`
- `scheme` optimizaiton stomatal model scheme OSMWang

This function uses the Wang et al. (2020) model. `∂Θ/∂E = A / (e_crit - e_leaf)`, unit `[mol mol⁻¹]`.

"""
function get_marginal_penalty(canopyi::CanopyLayer, indx::Number, scheme::OSMWang)
    FTYP = eltype(canopyi.an_list[indx])
    ∂Θ∂E = canopyi.an_list[indx] / (canopyi.ec_list[indx] - canopyi.e_list[indx]) * FTYP(1e-6)
    return ∂Θ∂E
end




"""
    get_marginal_penalty(canopyi, scheme::OSMWang)

List of marginal water penalty for all the leaves, given
- `canopyi` A [`CanopyLayer`](@ref) in [`Canopy`](@ref) in a [`Tree`](@ref)
- `scheme` optimizaiton stomatal model scheme OSMWang

This function uses the Wang et al. (2020) model. `∂Θ/∂E = A / (e_crit - e_leaf)`, unit `[mol mol⁻¹]`.

"""
function get_marginal_penalty(canopyi::CanopyLayer, scheme::OSMWang)
    FTYP      = eltype(canopyi.an_list[1])
    list_∂Θ∂E = canopyi.an_list ./ (canopyi.ec_list .- canopyi.e_list) .* FTYP(1e-6)
    list_∂Θ∂E[ list_∂Θ∂E .<= 0 ] .= FTYP(Inf)
    return list_∂Θ∂E
end




"""
    get_marginal_penalty(canopyi, indx, scheme::OSMWAP)

Marginal water penalty (`∂Θ/∂E`) for opening the stomata for a leaf in canopy layer, given
- `canopyi` A [`CanopyLayer`](@ref) in [`Canopy`](@ref) in a [`Tree`](@ref)
- `indx` indx'th [`Leaf`](@ref) in the `canopyi.leaf_list`
- `scheme` optimizaiton stomatal model scheme OSMWAP

This function uses the Anderegg et al. (2018) model. `∂Θ/∂E = (2*a*P + b*P) / K`, unit `[mol mol⁻¹]`.

"""
function get_marginal_penalty(canopyi::CanopyLayer, indx::Number, scheme::OSMWAP)
    leaf = canopyi[indx]
    ∂Θ∂E = -(2 * scheme.a + scheme.b) * leaf.p_element[end] / leaf.k_element[end] * leaf.k_sla * 10
    return ∂Θ∂E
end




"""
    get_marginal_penalty(canopyi, scheme::OSMWAP)

List of marginal water penalty for all the leaves, given
- `canopyi` A [`CanopyLayer`](@ref) in [`Canopy`](@ref) in a [`Tree`](@ref)
- `scheme` optimizaiton stomatal model scheme OSMWAP

This function uses the Anderegg et al. (2018) model. `∂Θ/∂E = (2*a*P + b*P) / K`, unit `[mol mol⁻¹]`.

"""
function get_marginal_penalty(canopyi::CanopyLayer, scheme::OSMWAP)
    list_p    = [leaf.p_element[end] for leaf in canopyi.leaf_list]
    list_k    = [leaf.k_element[end] for leaf in canopyi.leaf_list]
    list_∂Θ∂E = -(2 * scheme.a + scheme.b) .* list_p ./ list_k .* canopyi.leaf_list[1].k_sla .* 10
    return list_∂Θ∂E
end




"""
    get_marginal_penalty(canopyi, indx, scheme::OSMWAPMod)

Marginal water penalty (`∂Θ/∂E`) for opening the stomata for a leaf in canopy layer, given
- `canopyi` A [`CanopyLayer`](@ref) in [`Canopy`](@ref) in a [`Tree`](@ref)
- `indx` indx'th [`Leaf`](@ref) in the `canopyi.leaf_list`
- `scheme` optimizaiton stomatal model scheme OSMWAPMod

This function is modified from Anderegg et al. (2018) model. `∂Θ/∂E = a * A * P / K`, unit `[mol mol⁻¹]`.

"""
function get_marginal_penalty(canopyi::CanopyLayer, indx::Number, scheme::OSMWAPMod)
    leaf = canopyi[indx]
    ∂Θ∂E = -scheme.a * canopyi.an_list[indx] * leaf.p_element[end] / leaf.k_element[end] * leaf.k_sla * 10
    return ∂Θ∂E
end




"""
    get_marginal_penalty(canopyi, scheme::OSMWAPMod)

List of marginal water penalty for all the leaves, given
- `canopyi` A [`CanopyLayer`](@ref) in [`Canopy`](@ref) in a [`Tree`](@ref)
- `scheme` optimizaiton stomatal model scheme OSMWAPMod

This function is modified from Anderegg et al. (2018) model. `∂Θ/∂E = a * A * P / K`, unit `[mol mol⁻¹]`.

"""
function get_marginal_penalty(canopyi::CanopyLayer, scheme::OSMWAPMod)
    list_p    = [leaf.p_element[end] for leaf in canopyi.leaf_list]
    list_k    = [leaf.k_element[end] for leaf in canopyi.leaf_list]
    list_∂Θ∂E = -scheme.a .* canopyi.an_list .* list_p ./ list_k .* canopyi.leaf_list[1].k_sla .* 10
    return list_∂Θ∂E
end
