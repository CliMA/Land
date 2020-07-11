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
            tree::TreeSimple{FT},
            ini::FT = FT(0.5)
) where {FT<:AbstractFloat}
    # calculate maximal flow for whole tree, remember to times leaf area
    _kr = (tree.root).k_max;
    _ks = (tree.stem).k_max;
    _kl = (tree.leaf).k_sla / (tree.leaf).f_vis * (tree.leaf).area;
    _kt = 1 / (1 / _kr + 1 / _ks + 1/ _kl);
    _fh = -1 * (tree.leaf).p_crt * _kt;
    _fl = FT(0);
    _fx = min((_fh+_fl)/2, ini);
    _ms = NewtonBisectionMethod{FT}(_fl, _fh, _fx);
    _rt = SolutionTolerance{FT}(1e-3);
    @inline f(x) = xylem_p_from_flow(tree, x) - (tree.leaf).p_crt;
    _solut  = find_zero(f, _ms, _rt);

    return _solut
end




function tree_e_crit(
            tree::AbstractPlantHS{FT},
            ini::FT = FT(0.5)
) where {FT<:AbstractFloat}
    # calculate maximal flow
    _hs = tree.leaves[1];
    _fh = (-_hs.p_crt) * _hs.k_sla / _hs.f_vis * _hs.area;
    _fl = FT(0);
    _fx = min((_fh+_fl)/2, ini);
    _ms = NewtonBisectionMethod{FT}(_fl, _fh, _fx);
    _rt = SolutionTolerance{FT}(1e-3);
    @inline f(x) = xylem_p_from_flow(tree, x) - tree.leaves[1].p_crt;
    _solut  = find_zero(f, _ms, _rt);

    return _solut
end








###############################################################################
#
# Whole plant hydraulics --- Element and whole-tree conductance ratio
#
###############################################################################
"""
    plant_conductances!(tree::TreeSimple{FT}) where {FT<:AbstractFloat}

Update plant total hydraulic information in each element and whole plant, given
- `tree` [`TreeSimple`](@ref) type struct
"""
function plant_conductances!(
            tree::TreeSimple{FT}
) where {FT<:AbstractFloat}
    # Calculate the resistances
    r_m_root = 1 / (tree.root).k_max;
    r_m_stem = 1 / (tree.stem).k_max;
    r_m_leaf = 1 / ((tree.leaf).k_sla * (tree.leaf).area);
    r_m_tree = r_m_root + r_m_stem + r_m_leaf

    # Calculate the relative resistances
    r_root = FT(0);
    r_stem = FT(0);
    r_leaf = FT(0);
    for i in 1:10
        r_root += (tree.root).k_history[i];
        r_stem += (tree.stem).k_history[i];
        r_leaf += (tree.leaf).k_history[i];
    end
    r_root /= 10;
    r_stem /= 10;
    r_leaf /= 10;
    r_tree  = r_root * r_m_root + r_stem * r_m_stem + r_leaf * r_m_leaf;

    # Update the kr info
    tree.krs[1] = 1 / r_root;
    tree.krs[2] = 1 / r_stem;
    tree.krs[3] = 1 / r_leaf;
    tree.krs[4] = r_m_tree / r_tree;

    return nothing
end
