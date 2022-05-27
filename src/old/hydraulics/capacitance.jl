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
#=
function update_PVF!(
            tree::MonoGrassSPAC{FT},
            Δt::FT
) where {FT<:AbstractFloat}
    @unpack cache_k, cache_p, cache_q, roots = tree;
    leaves = tree.leaves;

    # 0. note that leaf flow rates need to be updated outside this function
    # 1. update leaf PVF for each layer
    f_sum::FT = 0;
    for leaf in leaves
        update_PVF!(leaf, Δt);
        f_sum += (leaf).q_in * (leaf).area;
    end

    # 2. update flow rate for roots, partition total flow rates into roots
    update_PVF!(roots, cache_k, cache_p, cache_q, f_sum, Δt);

    # 3. update pressure profiles
    xylem_pressure_profile!(tree; update=false);

    return nothing
end




function update_PVF!(
            tree::MonoPalmSPAC{FT},
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
        f_sum += (leaf).q_in * (leaf).area;
    end

    # 2. update PVF for trunk
    (trunk).q_out = f_sum;
    update_PVF!(trunk, Δt);

    # 3. update flow rate for roots, partition total flow rates into roots
    update_PVF!(roots, cache_k, cache_p, cache_q, (trunk).q_in, Δt);

    # 4. update pressure profiles
    xylem_pressure_profile!(tree; update=false);

    return nothing
end




function update_PVF!(
            tree::MonoTreeSPAC{FT},
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
        stem.q_out = (leaf).q_in * (leaf).area;
        update_PVF!(stem, Δt);
        f_sum += (stem).q_in;
    end

    # 2. update PVF for trunk
    (trunk).q_out = f_sum;
    update_PVF!(trunk, Δt);

    # 3. update flow rate for roots, partition total flow rates into roots
    update_PVF!(roots, cache_k, cache_p, cache_q, (trunk).q_in, Δt);

    # 4. update pressure profiles
    xylem_pressure_profile!(tree; update=false);

    return nothing
end
=#
