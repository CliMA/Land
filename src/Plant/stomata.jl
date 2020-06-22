#=
###############################################################################
#
# Calculate the optimization model gain
# These functions passed the FT test
# These functions are documented in the Plant page
#
###############################################################################
"""
    get_marginal_gain(canopyi::CanopyLayer, indx::Number, photo_para_set::AbstractPhotoModelParaSet)
    get_marginal_gain(canopyi::CanopyLayer, photo_para_set::AbstractPhotoModelParaSet)

Marginal water use efficiency `∂A/∂E` for a leaf in a canopy layer, given
- `canopyi` A [`CanopyLayer`](@ref) in a [`Tree`](@ref)
- `indx` indx'th [`Leaf`](@ref) in the `canopyi.leaf_list`
- `photo_para_set` An `AbstractPhotoModelParaSet` type parameter set
"""
function get_marginal_gain(canopyi::CanopyLayer, indx::Number, photo_para_set::AbstractPhotoModelParaSet)
    FTYP      = eltype(canopyi.f_view)
    leaf_e1   = canopyi.e_list[indx] + FTYP(1E-6)
    leaf_d    = (saturation_vapor_pressure(canopyi.t_list[indx]) - canopyi.p_H₂O) /  canopyi.p_atm
    leaf_glw1 = leaf_e1 / leaf_d
    leaf_glc1 = leaf_glw1 / FTYP(1.6) / ( 1 + canopyi.g_ias_c * leaf_glw1^(canopyi.g_ias_e) )
    _envir    = EnvironmentConditions{FTYP}(p_a = canopyi.p_a, p_atm = canopyi.p_atm, p_O₂ = canopyi.p_O₂)
    anagrpi   = an_ag_r_pi_from_glc(
                                    photo_para_set,
                                    leaf_glc1,
                                    canopyi.v_max,
                                    canopyi.j_max,
                                    canopyi.p_max,
                                    canopyi.t_list[indx],
                                    canopyi.par_list[indx],
                                    canopyi.r_25,
                                    _envir)
    leaf_a1   = anagrpi[1]

    # mol CO₂ mol⁻¹ H₂O, Δe = 1e-6 mol (H₂O) and μmol = 1e-6 mol (CO₂) cancel out
    return leaf_a1 - canopyi.an_list[indx]
end

function get_marginal_gain(canopyi::CanopyLayer, photo_para_set::AbstractPhotoModelParaSet)
    FTYP      = eltype(canopyi.f_view)
    leaf_e1   = canopyi.e_list .+ FTYP(1E-6)
    leaf_glw1 = leaf_e1 ./ canopyi.d_list
    leaf_glc1 = leaf_glw1 ./ FTYP(1.6) ./ ( 1 .+ canopyi.g_ias_c .* canopyi.glw_list .^ (canopyi.g_ias_e) )
    _envir    = EnvironmentConditions{FTYP}(p_a = canopyi.p_a, p_atm = canopyi.p_atm, p_O₂ = canopyi.p_O₂)
    anagrpi_l = an_ag_r_pi_from_glc(
                                    photo_para_set,
                                    leaf_glc1,
                                    canopyi.v_max,
                                    canopyi.j_max,
                                    canopyi.p_max,
                                    canopyi.t_list,
                                    canopyi.par_list,
                                    canopyi.r_25,
                                    _envir)
    leaf_a1   = anagrpi_l[1]
    # mol CO₂ mol⁻¹ H₂O, 1e-6 mol (H₂O) and μmol = 1e-6 mol (CO₂) cancel out
    return leaf_a1 .- canopyi.an_list
end








###############################################################################
#
# Calculate the optimization model penalty
# These functions passed the FT test
# These functions are documented in the Plant page
#
###############################################################################
"""
    get_marginal_penalty(canopyi::CanopyLayer, indx::Number, scheme::AbstractOptimizationStomatalModel)
    get_marginal_penalty(canopyi::CanopyLayer, scheme::AbstractOptimizationStomatalModel)

Marginal water penalty (`∂Θ/∂E`) for opening the stomata for a leaf in canopy layer, given
- `canopyi` A [`CanopyLayer`](@ref) in a [`Tree`](@ref)
- `indx` indx'th [`Leaf`](@ref) in the `canopyi.leaf_list`
- `scheme` An `AbstractOptimizationStomatalModel` type optimizaiton stomatal model scheme

`OSMEller` scheme uses the Eller et al. (2018) model. `∂Θ/∂E = A / K * dK/dE`, unit `[mol mol⁻¹]`.
To make the Eller model suitable for multiple layer canopy, K at the leaf xylem end is used. Also, to make the model capable of handling drought history, loss of conductance is evaluated assuming no drought legacy from a more negative `p_leaf`.

`OSMSperry` uses the Sperry et al. (2017) model. `∂Θ/∂E = A_max / K_max * dK/dE`.
To make the Sperry model suitable for multiple layer canopy, K at the leaf xylem end is used. Also, to make the model capable of handling drought history, loss of conductance is evaluated assuming no drought legacy from a more negative `p_leaf`.

`OSMWang` function uses the Wang et al. (2020) model. `∂Θ/∂E = A / (e_crit - e_leaf)`.

`OSMWAP` uses the Anderegg et al. (2018) model. `∂Θ/∂E = (2*a*P + b) / K`.

`OSMWAPMod` uses modified Anderegg et al. (2018) model. `∂Θ/∂E = a * A * P / K`.
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
    ∂Θ∂E = (k_0 - k_1) / de * canopyi.an_list[indx] / leaf.k_element[end] * leaf.k_sla * 10
    return ∂Θ∂E
end

function get_marginal_penalty(canopyi::CanopyLayer, scheme::OSMEller)
    FTYP      = eltype(canopyi.la)
    list_∂Θ∂E = canopyi.an_list .* 0
    for indx in 1:length(canopyi.leaf_list)
        list_∂Θ∂E[indx] = get_marginal_penalty(canopyi, indx, scheme)
    end
    list_∂Θ∂E[ list_∂Θ∂E .<= 0 ] .= FTYP(Inf)
    return list_∂Θ∂E
end

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

function get_marginal_penalty(canopyi::CanopyLayer, scheme::OSMSperry)
    FTYP      = eltype(canopyi.la)
    list_∂Θ∂E = canopyi.an_list .* 0
    for indx in 1:length(canopyi.leaf_list)
        list_∂Θ∂E[indx] = get_marginal_penalty(canopyi, indx, scheme)
    end
    list_∂Θ∂E[ list_∂Θ∂E .<= 0 ] .= FTYP(Inf)
    return list_∂Θ∂E
end

function get_marginal_penalty(canopyi::CanopyLayer, indx::Number, scheme::OSMWang)
    FTYP = eltype(canopyi.an_list[indx])
    ∂Θ∂E = canopyi.an_list[indx] / (canopyi.ec_list[indx] - canopyi.e_list[indx]) * FTYP(1e-6)
    return ∂Θ∂E
end

function get_marginal_penalty(canopyi::CanopyLayer, scheme::OSMWang)
    FTYP      = eltype(canopyi.an_list[1])
    list_∂Θ∂E = canopyi.an_list ./ (canopyi.ec_list .- canopyi.e_list) .* FTYP(1e-6)
    list_∂Θ∂E[ list_∂Θ∂E .<= 0 ] .= FTYP(Inf)
    return list_∂Θ∂E
end

function get_marginal_penalty(canopyi::CanopyLayer, indx::Number, scheme::OSMWAP)
    leaf = canopyi.leaf_list[indx]
    ∂Θ∂E = -(2 * scheme.a * leaf.p_element[end] + scheme.b) / leaf.k_element[end] * leaf.k_sla * 10
    return ∂Θ∂E
end

function get_marginal_penalty(canopyi::CanopyLayer, scheme::OSMWAP)
    list_p    = [leaf.p_element[end] for leaf in canopyi.leaf_list]
    list_k    = [leaf.k_element[end] for leaf in canopyi.leaf_list]
    list_∂Θ∂E = -(2 .* scheme.a .* list_p .+ scheme.b) ./ list_k .* canopyi.leaf_list[1].k_sla .* 10
    return list_∂Θ∂E
end

function get_marginal_penalty(canopyi::CanopyLayer, indx::Number, scheme::OSMWAPMod)
    leaf = canopyi.leaf_list[indx]
    ∂Θ∂E = -scheme.a * canopyi.an_list[indx] * leaf.p_element[end] / leaf.k_element[end] * leaf.k_sla * 10
    return ∂Θ∂E
end

function get_marginal_penalty(canopyi::CanopyLayer, scheme::OSMWAPMod)
    list_p    = [leaf.p_element[end] for leaf in canopyi.leaf_list]
    list_k    = [leaf.k_element[end] for leaf in canopyi.leaf_list]
    list_∂Θ∂E = -scheme.a .* canopyi.an_list .* list_p ./ list_k .* canopyi.leaf_list[1].k_sla .* 10
    return list_∂Θ∂E
end
=#







#=
###############################################################################
#
# Update tree with time
# These functions passed the FT test
# These functions are documented in the Plant page
#
###############################################################################
"""
    update_tree_with_time!(tree::Tree,
                           Δt::FT,
                           scheme::AbstractOptimizationStomatalModel;
                           updating::Bool=false)
    update_tree_with_time!(tree::Tree,
                           Δt::FT,
                           scheme::ESMGentine;
                           updating::Bool=false)
    update_tree_with_time!(tree::Tree,
                           Δt::FT,
                           scheme::AbstractEmpiricalStomatalModel;
                           updating::Bool=false)

Update the tree flux information and hydraulic parameters, given
- `tree`        One [`Tree`](@ref) type struct
- `Δt`          Non-steady state ``g_\\text{sw}`` time in `[s]`
- `scheme`      An `AbstractOptimizationStomatalModel` type optimization stomatal model scheme
- `updating`    Hydraulic update option. If `true`, plant hydraulic parameters will be updated for the [`Tree`](@ref)

The function updates the non-steady state stomatal conductance by
```math
\\Delta g_\\text{sw} = \\text{factor} \\cdot \\left( \\dfrac{∂A}{∂E} - \\dfrac{∂Θ}{∂E} \\right) \\cdot Δt
```
for optimization stomatal models. The ``\\dfrac{∂A}{∂E}`` is the same for every stomatal control model, but ``\\dfrac{∂Θ}{∂E}`` differs among models.

The function updates the stomatal conductance via a `Δglw = factor * (gs_ss - gs_nss)` for empirical stomatal models like Ball-Berry, Gentine, Leuning, and Medlyn models. The `gs_mod` stands for the gs calculated from model, and `gs_nss` stands for gs at non-steady state.
"""
function update_tree_with_time!(tree::Tree,
                                Δt::FT,
                                scheme::AbstractOptimizationStomatalModel;
                                updating::Bool=false) where {FT}
    # 0. unpack necessary structs
    @unpack branch_list,canopy_list,root_list = tree

    # 1. use the glw from last time instant and update e and a first, because Tleaf was from the last instant
    for canopyi in canopy_list
        # 1.1 update the e and a for each leaf, per leaf area
        # TODO make sure canopyi.p_H₂O is lower than saturation_vapor_pressure(canopyi.t_air)
        N               = length(canopyi.glc_list)
        canopyi.d_list  = (saturation_vapor_pressure(canopyi.t_list) .- canopyi.p_H₂O) ./ canopyi.p_atm
        canopyi.e_list  = canopyi.glw_list .* canopyi.d_list
        canopyi.q_list  = canopyi.e_list   .* canopyi.la_list
        _envir          = EnvironmentConditions{FT}(p_a = canopyi.p_a, p_atm = canopyi.p_atm, p_O₂ = canopyi.p_O₂)
        anagrpi_lists   = an_ag_r_pi_from_glc(
                                              tree.photo_para_set,
                                              canopyi.glc_list,
                                              canopyi.v_max,
                                              canopyi.j_max,
                                              canopyi.p_max,
                                              canopyi.t_list,
                                              canopyi.par_list,
                                              canopyi.r_25,
                                              _envir)
        canopyi.an_list = anagrpi_lists[1]
        canopyi.ag_list = anagrpi_lists[2]
        canopyi.r_list  = anagrpi_lists[3]
        canopyi.pi_list = anagrpi_lists[4]
    end

    # update the pressure profile in the trunk, branch, and leaf
    if updating
        # 2. update the pressure profiles for roots
        q_canopy_list      = [sum(canopyi.q_list) for canopyi in canopy_list]
        q_sum              = sum(q_canopy_list)
        p_base,q_root_list = get_p_base_q_list_from_q(tree, q_sum)
        for indx in 1:length(root_list)
            rooti = root_list[indx]
            update_struct_from_q!(rooti, q_root_list[indx])
        end

        # 3. update the pressure profile in the trunk
        tree.trunk.p_ups = p_base
        update_struct_from_q!(tree.trunk, q_sum)

        # 4. update the pressure profiles in the branches and leaves
        for indx in 1:length(branch_list)
            # 4.1 update the pressure profile in the branch
            branchi       = branch_list[indx]
            canopyi       = canopy_list[indx]
            branchi.p_ups = tree.trunk.p_dos
            update_struct_from_q!(branchi, q_canopy_list[indx])
            
            # 4.2 update the pressure profiles for leaves
            for ind_leaf in 1:length(canopyi.leaf_list)
                leaf = canopyi.leaf_list[ind_leaf]
                leaf.p_ups = branchi.p_dos
                update_struct_from_q!(leaf, canopyi.e_list[ind_leaf])
            end
        end
    end

    # 5. determine how much glw and glc should change with time
    for canopyi in canopy_list
        # 5.1 compute the ∂A∂E and ∂Θ∂E for each leaf
        list_∂A∂E = get_marginal_gain(canopyi, tree.photo_para_set)
        list_∂Θ∂E = get_marginal_penalty(canopyi, scheme)

        # update the glw and glc for each leaf
        glw_list = canopyi.glw_list + canopyi.gs_nssf .* (list_∂A∂E - list_∂Θ∂E) .* Δt
        glw_min  = canopyi.g_min .* relative_diffusive_coefficient(canopyi.t_list)
        glw_list = max.( glw_list, glw_min )

        canopyi.glw_list = glw_list
        canopyi.glc_list = canopyi.glw_list ./ FT(1.6) ./ ( 1 .+ canopyi.g_ias_c .* canopyi.glw_list .^ (canopyi.g_ias_e) )
    end
end

function update_tree_with_time!(tree::Tree, Δt::FT, scheme::ESMGentine; updating::Bool=false) where {FT}
    # 0. unpack necessary structs
    @unpack branch_list,canopy_list,root_list   = tree

    # 1. use the glw from last time instant and update e and a first, because Tleaf was from the last instant
    for canopyi in canopy_list
        # 1.1 update the e and a for each leaf, per leaf area
        # TODO make sure canopyi.p_H₂O is lower than saturation_vapor_pressure(canopyi.t_air)
        N               = length(canopyi.glc_list)
        canopyi.d_list  = (saturation_vapor_pressure(canopyi.t_list) .- canopyi.p_H₂O) ./ canopyi.p_atm
        canopyi.e_list  = canopyi.glw_list .* canopyi.d_list
        canopyi.q_list  = canopyi.e_list   .* canopyi.la_list
        _envir          = EnvironmentConditions{FT}(p_a = canopyi.p_a, p_atm = canopyi.p_atm, p_O₂ = canopyi.p_O₂)
        anagrpi_lists   = an_ag_r_pi_from_glc(
                                              tree.photo_para_set,
                                              canopyi.glc_list,
                                              canopyi.v_max,
                                              canopyi.j_max,
                                              canopyi.p_max,
                                              canopyi.t_list,
                                              canopyi.par_list,
                                              canopyi.r_25,
                                              _envir)
        canopyi.an_list = anagrpi_lists[1]
        canopyi.ag_list = anagrpi_lists[2]
        canopyi.r_list  = anagrpi_lists[3]
        canopyi.pi_list = anagrpi_lists[4]
    end

    # update the pressure profile in the trunk, branch, and leaf
    if updating
        # 2. update the pressure profiles for roots
        q_canopy_list      = [sum(canopyi.q_list) for canopyi in canopy_list]
        q_sum              = sum(q_canopy_list)
        p_base,q_root_list = get_p_base_q_list_from_q(tree, q_sum)
        for indx in 1:length(root_list)
            rooti = root_list[indx]
            update_struct_from_q!(rooti, q_root_list[indx])
        end

        # 3. update the pressure profile in the trunk
        tree.trunk.p_ups = p_base
        update_struct_from_q!(tree.trunk, q_sum)

        # 4. update the pressure profiles in the branches and leaves
        for indx in 1:length(branch_list)
            # 4.1 update the pressure profile in the branch
            branchi       = branch_list[indx]
            canopyi       = canopy_list[indx]
            branchi.p_ups = tree.trunk.p_dos
            update_struct_from_q!(branchi, q_canopy_list[indx])
            
            # 4.2 update the pressure profiles for leaves
            for ind_leaf in 1:length(canopyi.leaf_list)
                leaf = canopyi.leaf_list[ind_leaf]
                leaf.p_ups = branchi.p_dos
                update_struct_from_q!(leaf, canopyi.e_list[ind_leaf])
            end
        end
    end

    # 5. determine how much glw and glc should change with time
    for canopyi in canopy_list
        N = length(canopyi.glc_list)
        # 5.1 compute the gs at steady state, function to be added
        # TODO add gb and gm here for empirical model
        # Distinguidh p_s and p_a later
        k_ratio  = weibull_k_leaf_ratio(canopyi.leaf_list, canopyi.t_list)
        gss_list = get_empirical_gsw_from_model(scheme, canopyi.an_list, canopyi.p_atm, canopyi.p_a, k_ratio)

        # update the glw and glc for each leaf
        glw_list = canopyi.glw_list + canopyi.gs_nssf .* (gss_list - canopyi.glw_list) .* Δt
        glw_min  = canopyi.g_min .* relative_diffusive_coefficient(canopyi.t_list)
        glw_list = max.( glw_list, glw_min )

        canopyi.glw_list = glw_list
        canopyi.glc_list = canopyi.glw_list ./ FT(1.6) ./ ( 1 .+ canopyi.g_ias_c .* canopyi.glw_list .^ (canopyi.g_ias_e) )
    end
end

function update_tree_with_time!(tree::Tree, Δt::FT, scheme::AbstractEmpiricalStomatalModel; updating::Bool=false) where {FT}
    # 0. unpack necessary structs
    @unpack branch_list,canopy_list,root_list = tree

    # 1. use the glw from last time instant and update e and a first, because Tleaf was from the last instant
    for canopyi in canopy_list
        # 1.1 update the e and a for each leaf, per leaf area
        # TODO make sure canopyi.p_H₂O is lower than saturation_vapor_pressure(canopyi.t_air)
        N               = length(canopyi.glc_list)
        canopyi.d_list  = (saturation_vapor_pressure(canopyi.t_list) .- canopyi.p_H₂O) ./ canopyi.p_atm
        canopyi.e_list  = canopyi.glw_list .* canopyi.d_list
        canopyi.q_list  = canopyi.e_list   .* canopyi.la_list
        _envir          = EnvironmentConditions{FT}(p_a = canopyi.p_a, p_atm = canopyi.p_atm, p_O₂ = canopyi.p_O₂)
        anagrpi_lists   = an_ag_r_pi_from_glc(
                                              tree.photo_para_set,
                                              canopyi.glc_list,
                                              canopyi.v_max,
                                              canopyi.j_max,
                                              canopyi.p_max,
                                              canopyi.t_list,
                                              canopyi.par_list,
                                              canopyi.r_25,
                                              _envir)
        canopyi.an_list = anagrpi_lists[1]
        canopyi.ag_list = anagrpi_lists[2]
        canopyi.r_list  = anagrpi_lists[3]
        canopyi.pi_list = anagrpi_lists[4]
    end

    # update the pressure profile in the trunk, branch, and leaf
    if updating
        # 2. update the pressure profiles for roots
        q_canopy_list      = [sum(canopyi.q_list) for canopyi in canopy_list]
        q_sum              = sum(q_canopy_list)
        p_base,q_root_list = get_p_base_q_list_from_q(tree, q_sum)
        for indx in 1:length(root_list)
            rooti = root_list[indx]
            update_struct_from_q!(rooti, q_root_list[indx])
        end

        # 3. update the pressure profile in the trunk
        tree.trunk.p_ups = p_base
        update_struct_from_q!(tree.trunk, q_sum)

        # 4. update the pressure profiles in the branches and leaves
        for indx in 1:length(branch_list)
            # 4.1 update the pressure profile in the branch
            branchi       = branch_list[indx]
            canopyi       = canopy_list[indx]
            branchi.p_ups = tree.trunk.p_dos
            update_struct_from_q!(branchi, q_canopy_list[indx])
            
            # 4.2 update the pressure profiles for leaves
            for ind_leaf in 1:length(canopyi.leaf_list)
                leaf = canopyi.leaf_list[ind_leaf]
                leaf.p_ups = branchi.p_dos
                update_struct_from_q!(leaf, canopyi.e_list[ind_leaf])
            end
        end
    end

    # 5. determine how much glw and glc should change with time
    for canopyi in canopy_list
        N = length(canopyi.glc_list)
        # 5.1 compute the gs at steady state, function to be added
        if tree.photo_type=="C3"
            Γ_star = photo_TD_from_set(tree.photo_para_set.ΓsT, canopyi.t_list)
        else
            Γ_star = canopyi.t_list .* 0
        end
        rh   = canopyi.p_H₂O / saturation_vapor_pressure(canopyi.t_air)
        vpd  = saturation_vapor_pressure(canopyi.t_list) .- canopyi.p_H₂O
        beta = ones(FT, length(canopyi.t_air))
        # TODO add beta function later
        # TODO add gb and gm as well here
        gss_list = get_empirical_gsw_from_model(scheme, canopyi.an_list, canopyi.p_atm, canopyi.p_a, rh, canopyi.d_list, Γ_star, beta)

        # update the glw and glc for each leaf
        glw_list = canopyi.glw_list + canopyi.gs_nssf .* (gss_list - canopyi.glw_list) .* Δt
        glw_min  = canopyi.g_min .* relative_diffusive_coefficient(canopyi.t_list)
        glw_list = max.( glw_list, glw_min )

        canopyi.glw_list = glw_list
        canopyi.glc_list = canopyi.glw_list ./ FT(1.6) ./ ( 1 .+ canopyi.g_ias_c .* canopyi.glw_list .^ (canopyi.g_ias_e) )
    end
end
=#