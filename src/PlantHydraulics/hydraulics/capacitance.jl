###############################################################################
#
# Update xylem pressure and capacitance
#
###############################################################################
"""
    update_PVF!(hs::LeafHydraulics{FT}, Δt::FT) where {FT<:AbstractFloat}
    update_PVF!(hs::StemHydraulics{FT}, Δt::FT) where {FT<:AbstractFloat}
    update_PVF!(hs::RootHydraulics{FT}, Δt::FT) where {FT<:AbstractFloat}
    update_PVF!(hs::RootHydraulics{FT},
                Δt::FT,
                nss::Bool
    ) where {FT<:AbstractFloat}
    update_PVF!(roots::Array{RootHydraulics{FT},1},
                ks::Array{FT,1},
                ps::Array{FT,1},
                qs::Array{FT,1},
                q_sum::FT,
                Δt::FT
    ) where {FT<:AbstractFloat}
    update_PVF!(tree::GrassLikeOrganism{FT}, Δt::FT) where {FT<:AbstractFloat}
    update_PVF!(tree::PalmLikeOrganism{FT}, Δt::FT) where {FT<:AbstractFloat}
    update_PVF!(tree::TreeLikeOrganism{FT}, Δt::FT) where {FT<:AbstractFloat}

Update pressure, capacitance, and flow rates in hydraulic system, given
- `hs` [`AbstractHydraulicOrgan`](@ref) or [`AbstractPlantOrganism`](@ref) type
    struct
- `Δt` Time interval
- `nss` NonSteadyStateMode placeholder (useless)
- `ks` Container for conductance in each root layer
- `ps` Container for end xylem pressure in each layer
- `qs` Container for flow rate out of each layer
- `q_sum` Total flow rate out of the roots
- `tree` [`GrassLikeOrganism`](@ref), [`PalmLikeOrganism`](@ref), or
    [`TreeLikeOrganism`](@ref) type struct

Note that this function only updates the equilibrium pressure in the tissue,
    but not the xylem flow pressure. The difference between the two pressures
    is used to drive water exchange between xylem and capacictance tissues.
"""
function update_PVF!(hs::LeafHydraulics{FT}, Δt::FT) where {FT<:AbstractFloat}
    # unpack values
    @unpack f_vis, p_element, p_leaf, q_out, pv, T_sap, v_maximum,
            v_storage = hs;

    # calculate the flow rate out of the capacitance
    # use buffer_rate here to avoid alloaction issues
    f_cap = (hs.p_storage - p_leaf) * buffer_rate(pv) / f_vis * v_maximum;
    if (f_cap > 0) && (hs.v_storage <= f_cap * Δt)
        f_cap = hs.v_storage / Δt;
    end

    # update flow into the tissue, storage and the tissue pressure (p_storage)
    hs.q_in       = q_out - f_cap;
    hs.v_storage -= f_cap * Δt;
    hs.p_storage  = p_from_volume(pv, hs.v_storage/v_maximum, T_sap);

    return nothing
end




function update_PVF!(hs::StemHydraulics{FT}, Δt::FT) where {FT<:AbstractFloat}
    # unpack values
    @unpack f_vis, N, q_element, p_element, p_storage, pv, q_out, T_sap,
            v_maximum, v_storage = hs;

    # update flow rate in, storage volume and pressure per slice
    f_sum::FT = 0;
    for i in N:-1:1
        f_cap = (p_storage[i] - p_element[i]) * buffer_rate(pv) / f_vis *
                v_maximum[i];
        if (f_cap > 0) && (v_storage[i] <= f_cap * Δt)
            f_cap = v_storage[i] / Δt;
        end
        v_storage[i] -= f_cap * Δt;
        p_storage[i]  = p_from_volume(pv, v_storage[i]/v_maximum[i], T_sap);
        q_element[i]  = q_out - f_sum;
        f_sum        += f_cap;
    end
    hs.q_in = q_out - f_sum;

    return nothing
end




function update_PVF!(hs::RootHydraulics{FT}, Δt::FT) where {FT<:AbstractFloat}
    # unpack values
    @unpack f_vis, N, p_element, p_storage, pv, q_buffer, q_diff, q_element,
            q_out, T_sap, v_maximum, v_storage = hs;

    # update flow rate in, storage volume and pressure per slice
    f_sum::FT = 0;
    for i in N:-1:1
        f_cap = (p_storage[i] - p_element[i]) * buffer_rate(pv) / f_vis *
                v_maximum[i];
        if (f_cap > 0) && (v_storage[i] <= f_cap * Δt)
            f_cap = v_storage[i] / Δt;
        end
        q_buffer[i]   = f_cap;
        q_diff[i]     = f_sum;
        v_storage[i] -= f_cap * Δt;
        p_storage[i]  = p_from_volume(pv, v_storage[i]/v_maximum[i], T_sap);
        q_element[i]  = q_out - f_sum;
        f_sum        += f_cap;
    end
    hs.q_in = q_out - f_sum;

    return nothing
end




function update_PVF!(
            hs::RootHydraulics{FT},
            Δt::FT,
            nss::Bool
) where {FT<:AbstractFloat}
    # unpack values
    @unpack f_vis, N, p_element, p_storage, pv, q_buffer, q_diff, q_element,
            q_out, v_maximum, v_storage = hs;

    # update flow rate in, storage volume and pressure per slice
    f_sum::FT = 0;
    for i in N:-1:1
        f_cap = (p_storage[i] - p_element[i]) * buffer_rate(pv) /f_vis *
                v_maximum[i];
        if (f_cap > 0) && (v_storage[i] <= f_cap * Δt)
            f_cap = v_storage[i] / Δt;
        end
        q_buffer[i] = f_cap;
        q_diff[i]   = f_sum;
        f_sum      += f_cap;
    end

    return nothing
end




function update_PVF!(
            roots::Array{RootHydraulics{FT},1},
            cache_k::Array{FT,1},
            cache_p::Array{FT,1},
            cache_q::Array{FT,1},
            q_sum::FT,
            Δt::FT
) where {FT<:AbstractFloat}
    # update flow rate for roots, partition total flow rates into roots
    # algorithm from function `roots_flow!`
    # 1. update root buffer rate
    for root in roots
        update_PVF!(root, Δt, true);
    end

    # 2. calculate the flow rate in each root layer
    count = 0;
    while true
        count += 1;
        if count > 20
            break
        end

        # 2.1 flush the values to ks, ps, and qs
        for i in eachindex(roots)
            root  = roots[i];
            _q_in = cache_q[i] - sum(root.q_buffer);
            root.q_element .= cache_q[i] .- root.q_diff;
            _p,_k = root_pk(root, _q_in, root.q_element);
            cache_p[i] = _p;
            cache_k[i] = _k;
        end

        # 2.2 use ps and ks to compute the Δq to adjust
        pm = mean(cache_p);
        for i in eachindex(roots)
            cache_q[i] -= (pm - cache_p[i]) * cache_k[i]
        end

        # 2.3 adjust the qs so that sum(qs) = f_sum
        q_diff = sum(cache_q) - q_sum;
        if (abs(q_diff) < max(FT(1e-6), eps(FT))) &&
           (maximum(cache_p) - minimum(cache_p) < 1e-4)
            break
        end
        k_sum  = sum(cache_k);
        for i in eachindex(roots)
            cache_q[i] -= q_diff * cache_k[1] / k_sum;
        end
    end

    # 2.4 update root PVF again
    for i in eachindex(roots)
        root = roots[i];
        root.q_out = cache_q[i];
        update_PVF!(root, Δt);
    end

    return nothing
end




function update_PVF!(
            tree::GrassLikeOrganism{FT},
            Δt::FT
) where {FT<:AbstractFloat}
    @unpack cache_k, cache_p, cache_q, roots = tree;
    leaves = tree.leaves;

    # 0. note that leaf flow rates need to be updated outside this function
    # 1. update leaf PVF for each layer
    f_sum::FT = 0;
    for leaf in leaves
        update_PVF!(leaf, Δt);
        f_sum += leaf.q_in * leaf.area;
    end

    # 2. update flow rate for roots, partition total flow rates into roots
    update_PVF!(roots, cache_k, cache_p, cache_q, f_sum, Δt);

    # 3. update pressure profiles
    pressure_profile!(tree, NonSteadyStateMode(); update=false);

    return nothing
end




function update_PVF!(
            tree::PalmLikeOrganism{FT},
            Δt::FT
) where {FT<:AbstractFloat}
    @unpack cache_k, cache_p, cache_q, roots = tree;
    leaves = tree.leaves;
    trunk = tree.trunk;

    # 0. note that leaf flow rates need to be updated outside this function
    # 1. update leaf PVF for each layer
    f_sum::FT = 0;
    for leaf in leaves
        update_PVF!(leaf, Δt);
        f_sum += leaf.q_in * leaf.area;
    end

    # 2. update PVF for trunk
    trunk.q_out = f_sum;
    update_PVF!(trunk, Δt);

    # 3. update flow rate for roots, partition total flow rates into roots
    update_PVF!(roots, cache_k, cache_p, cache_q, trunk.q_in, Δt);

    # 4. update pressure profiles
    pressure_profile!(tree, NonSteadyStateMode(); update=false);

    return nothing
end




function update_PVF!(
            tree::TreeLikeOrganism{FT},
            Δt::FT
) where {FT<:AbstractFloat}
    @unpack branch, cache_k, cache_p, cache_q, leaves, roots = tree;
    trunk = tree.trunk;

    # 0. note that leaf flow rates need to be updated outside this function
    # 1. update leaf and stem PVF for each layer
    f_sum::FT = 0;
    for i in eachindex(leaves)
        leaf = leaves[i];
        stem = branch[i];
        update_PVF!(leaf, Δt);
        stem.q_out = leaf.q_in * leaf.area;
        update_PVF!(stem, Δt);
        f_sum += stem.q_in;
    end

    # 2. update PVF for trunk
    trunk.q_out = f_sum;
    update_PVF!(trunk, Δt);

    # 3. update flow rate for roots, partition total flow rates into roots
    update_PVF!(roots, cache_k, cache_p, cache_q, trunk.q_in, Δt);

    # 4. update pressure profiles
    pressure_profile!(tree, NonSteadyStateMode(); update=false);

    return nothing
end
