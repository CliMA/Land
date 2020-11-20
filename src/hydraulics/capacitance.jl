###############################################################################
#
# Capacitance buffer flow
#
###############################################################################
"""
    buffer_rate(pv::AbstractCapacity{FT}) where {FT<:AbstractFloat}

Return the buffer rate, given
- `pv` [`AbstractCapacity`](@ref type struct)
"""
function buffer_rate(
            pv::PVCurveLinear{FT}
) where {FT<:AbstractFloat}
    return pv.k_refill
end








###############################################################################
#
# Update xylem pressure and capacitance
#
###############################################################################
"""
    update_PVF!(hs::AbstractHydraulicSystem{FT},
                Δt::FT
    ) where {FT<:AbstractFloat}
    update_PVF!(hs::AbstractPlantHS{FT},
                Δt::FT
    ) where {FT<:AbstractFloat}

Update pressure, capacitance, and flow rates in hydraulic system, given
- `hs` [`AbstractHydraulicSystem`](@ref) or [`AbstractPlantHS`](@ref) type
    struct
- `Δt` Time interval

Note that this function only updates the equilibrium pressure in the tissue,
    but not the xylem flow pressure. The difference between the two pressures
    is used to drive water exchange between xylem and capacictance tissues.
"""
function update_PVF!(hs::LeafHydraulics{FT}, Δt::FT) where {FT<:AbstractFloat}
    # unpack values
    @unpack p_element, p_leaf, q_out, pv, v_maximum, v_storage = hs;

    # calculate the flow rate out of the capacitance
    # use buffer_rate here to avoid alloaction issues
    f_cap = (hs.p_storage - p_leaf) * buffer_rate(pv) * v_maximum;
    if (f_cap > 0) && (hs.v_storage <= f_cap * Δt)
        f_cap = hs.v_storage / Δt;
    end

    # update flow into the tissue, storage and the tissue pressure (p_storage)
    hs.q_in       = q_out - f_cap;
    hs.v_storage -= f_cap * Δt;
    hs.p_storage  = p_from_volume(pv, hs.v_storage/v_maximum);

    return nothing
end




function update_PVF!(hs::StemHydraulics{FT}, Δt::FT) where {FT<:AbstractFloat}
    # unpack values
    @unpack q_element, p_element, p_storage, pv, q_out, v_maximum, v_storage =
            hs;

    # update flow rate in, storage volume and pressure per slice
    f_sum::FT = 0;
    for i in 10:-1:1
        f_cap = (p_storage[i] - p_element[i]) * buffer_rate(pv) * v_maximum[i];
        if (f_cap > 0) && (v_storage[i] <= f_cap * Δt)
            f_cap = v_storage[i] / Δt;
        end
        v_storage[i] -= f_cap * Δt;
        p_storage[i]  = p_from_volume(pv, v_storage[i]/v_maximum[i]);
        q_element[i]  = q_out - f_sum;
        f_sum        += f_cap;
    end
    hs.q_in = q_out - f_sum;

    return nothing
end




function update_PVF!(
            hs::RootHydraulics{FT},
            Δt::FT,
            nss::Bool = true
) where {FT<:AbstractFloat}
    # unpack values
    @unpack p_element, p_storage, pv, q_buffer, q_diff, q_element, q_out,
            v_maximum, v_storage = hs;

    # update flow rate in, storage volume and pressure per slice
    f_sum::FT = 0;
    for i in 10:-1:1
        f_cap = (p_storage[i] - p_element[i]) * buffer_rate(pv) * v_maximum[i];
        if (f_cap > 0) && (v_storage[i] <= f_cap * Δt)
            f_cap = v_storage[i] / Δt;
        end
        q_buffer[i] = f_cap;
        q_diff[i]   = f_sum;
        # update these only when nss is false
        if !nss
            v_storage[i] -= f_cap * Δt;
            p_storage[i]  = p_from_volume(pv, v_storage[i]/v_maximum[i]);
            q_element[i]  = q_out - f_sum;
        end
        f_sum        += f_cap;
    end
    # update q_in only when nss is false
    if !nss
        hs.q_in = q_out - f_sum;
    end

    return nothing
end




function update_PVF!(
            roots::Array{RootHydraulics{FT},1},
            container_k::Array{FT,1},
            container_p::Array{FT,1},
            container_q::Array{FT,1},
            q_sum::FT,
            Δt::FT
) where {FT<:AbstractFloat}
    # update flow rate for roots, partition total flow rates into roots
    # algorithm from function `recalculate_roots_flow!`
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
            _q_in = container_q[i] - sum(root.q_buffer);
            root.q_element .= container_q[i] .- root.q_diff;
            _p,_k = root_pk_from_flow(root, _q_in, root.q_element);
            container_p[i] = _p;
            container_k[i] = _k;
        end

        # 2.2 use ps and ks to compute the Δq to adjust
        pm = mean(container_p);
        for i in eachindex(roots)
            container_q[i] -= (pm - container_p[i]) * container_k[i]
        end

        # 2.3 adjust the qs so that sum(qs) = f_sum
        q_diff = sum(container_q) - q_sum;
        if abs(q_diff) < max(FT(1e-6), eps(FT))
            break
        end
        k_sum  = sum(container_k);
        for i in eachindex(roots)
            container_q[i] -= q_diff * container_k[1] / k_sum;
        end
    end

    # 2.4 update root PVF again
    for i in eachindex(roots)
        root = roots[i];
        root.q_out = container_q[i];
        update_PVF!(root, Δt, false);
    end

    return nothing
end




function update_PVF!(tree::GrassLikeHS{FT}, Δt::FT) where {FT<:AbstractFloat}
    @unpack container_k, container_p, container_q, leaves, roots = tree;

    # 0. note that leaf flow rates need to be updated outside this function
    # 1. update leaf PVF for each layer
    f_sum::FT = 0;
    for leaf in leaves
        update_PVF!(leaf, Δt);
        f_sum += leaf.q_in * leaf.area;
    end

    # 2. update flow rate for roots, partition total flow rates into roots
    update_PVF!(roots, container_k, container_p, container_q, f_sum, Δt);

    # 3. update pressure profiles
    pressure_profile!(tree; update=false);

    return nothing
end




function update_PVF!(tree::PalmLikeHS{FT}, Δt::FT) where {FT<:AbstractFloat}
    @unpack container_k, container_p, container_q, leaves, roots, trunk = tree;

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
    update_PVF!(roots, container_k, container_p, container_q, trunk.q_in, Δt);

    # 4. update pressure profiles
    pressure_profile!(tree; update=false);

    return nothing
end




function update_PVF!(tree::TreeLikeHS{FT}, Δt::FT) where {FT<:AbstractFloat}
    @unpack branch, container_k, container_p, container_q, leaves, roots,
            trunk = tree;

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
    update_PVF!(roots, container_k, container_p, container_q, trunk.q_in, Δt);

    # 4. update pressure profiles
    pressure_profile!(tree; update=false);

    return nothing
end
