###############################################################################
#
# Whole plant hydraulics --- Tree end pressure from flow
#
###############################################################################
"""
    tree_p_from_flow(tree::GrassLikeHS{FT}, flow::FT) where {FT<:AbstractFloat}
    tree_p_from_flow(tree:: PalmLikeHS{FT}, flow::FT) where {FT<:AbstractFloat}
    tree_p_from_flow(tree:: TreeLikeHS{FT}, flow::FT) where {FT<:AbstractFloat}

Calculate the leaf xylem end pressure, given
- `tree` [`GrassLikeHS`](@ref) or [`PalmLikeHS`](@ref) or [`TreeLikeHS`](@ref)
    type hydraulics system
- `flow` Given flow rate
"""
function tree_p_from_flow(
            tree::GrassLikeHS{FT},
            flow::FT
) where {FT<:AbstractFloat}
    # Total leaf area
    tla = sum([tree.leaves[i].area for i in 1:tree.n_canopy])

    if tree.n_root == tree.n_canopy == 1
        # calculate the p_dos for roots
        p_dos = xylem_p_from_flow(tree.roots[1], flow);
        tree.leaves[1].p_ups = p_dos;

        # calculate the p_dos for leaves
        p_dos = xylem_p_from_flow(tree.leaves[1], flow/tla);

        return p_dos
    elseif tree.n_canopy == 1
        # calculate the p_dos for roots
        _qs,_p = root_qs_p_from_q(tree.roots, flow, (tree.leaves[1]).p_ups);
        tree.leaves[1].p_ups = _p;

        # calculate the p_dos for leaves
        p_dos = xylem_p_from_flow(tree.leaves[1], flow/tla);

        return p_dos
    else
        println("No function applicable for multi-layer canopy!")
        return ErrorException("Error!");
    end
end




function tree_p_from_flow(
            tree::PalmLikeHS{FT},
            flow::FT
) where {FT<:AbstractFloat}
    # Total leaf area
    tla = sum([tree.leaves[i].area for i in 1:tree.n_canopy])

    if tree.n_root == tree.n_canopy == 1
        # calculate the p_dos for roots
        p_dos = xylem_p_from_flow(tree.roots[1], flow);
        (tree.trunk).p_ups = p_dos;

        # calculate the p_dos for trunk
        p_dos = xylem_p_from_flow(tree.trunk, flow);
        tree.leaves[1].p_ups = p_dos;

        # calculate the p_dos for leaves
        p_dos = xylem_p_from_flow(tree.leaves[1], flow/tla);

        return p_dos
    elseif tree.n_canopy == 1
        # calculate the p_dos for roots
        _qs,_p = root_qs_p_from_q(tree.roots, flow, (tree.trunk).p_ups);
        (tree.trunk).p_ups = _p;

        # calculate the p_dos for trunk
        p_dos = xylem_p_from_flow(tree.trunk, flow);
        tree.leaves[1].p_ups = p_dos;

        # calculate the p_dos for leaves
        p_dos = xylem_p_from_flow(tree.leaves[1], flow/tla);

        return p_dos
    else
        println("No function applicable for multi-layer canopy!")
        return ErrorException("Error!");
    end
end




function tree_p_from_flow(
            tree::TreeLikeHS{FT},
            flow::FT
) where {FT<:AbstractFloat}
    # Total leaf area
    tla = sum([tree.leaves[i].area for i in 1:tree.n_canopy])

    if tree.n_root == tree.n_canopy == 1
        # calculate the p_dos for roots
        p_dos = xylem_p_from_flow(tree.roots[1], flow);
        (tree.trunk).p_ups = p_dos;

        # calculate the p_dos for trunk
        p_dos = xylem_p_from_flow(tree.trunk, flow);
        tree.branch[1].p_ups = p_dos;

        # calculate the p_dos for branch
        p_dos = xylem_p_from_flow(tree.branch[1], flow);
        tree.leaves[1].p_ups = p_dos;

        # calculate the p_dos for leaves
        p_dos = xylem_p_from_flow(tree.leaves[1], flow/tla);

        return p_dos
    elseif tree.n_canopy == 1
        # calculate the p_dos for roots
        _qs,_p = root_qs_p_from_q(tree.roots, flow, (tree.trunk).p_ups);
        (tree.trunk).p_ups = _p;

        # calculate the p_dos for trunk
        p_dos = xylem_p_from_flow(tree.trunk, flow);
        tree.branch[1].p_ups = p_dos;

        # calculate the p_dos for branch
        p_dos = xylem_p_from_flow(tree.branch[1], flow);
        tree.leaves[1].p_ups = p_dos;

        # calculate the p_dos for leaves
        p_dos = xylem_p_from_flow(tree.leaves[1], flow/tla);

        return p_dos
    else
        println("No function applicable for multi-layer canopy!")
        return ErrorException("Error!");
    end
end








###############################################################################
#
# Whole plant hydraulics --- Tree critical pressure
#
###############################################################################
"""
    tree_e_crit(tree::PalmLikeHS{FT}, ini::FT) where {FT<:AbstractFloat}

Calculate the critical flow rate for the tree per LA, given
- `tree` [`GrassLikeHS`](@ref) or [`PalmLikeHS`](@ref) or [`TreeLikeHS`](@ref)
    type hydraulics system
- `ini` Initial guess of the critical flow rate
"""
function tree_e_crit(
            tree::AbstractPlantHS{FT},
            ini::FT = FT(0.5)
) where {FT<:AbstractFloat}
    # create a tree copy to avoid changing the tree
    # tree_copy = deepcopy(tree);

    # calculate maximal flow
    _hs = tree.leaves[1];
    _fh = (-_hs.p_crt) * _hs.k_sla / _hs.f_vis;
    _fl = FT(0);
    _fx = min((_fh+_fl)/2, ini);
    _ms = NewtonBisectionMethod{FT}(_fl, _fh, _fx);
    _rt = ResidualTolerance{FT}(1e-5);
    @inline f(x) = tree_p_from_flow(tree, x) - tree.leaves[1].p_crt;
    _solut  = find_zero(f, _ms, _rt);
    _ec = _solut;

    return _ec
end
