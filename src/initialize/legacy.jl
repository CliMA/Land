###############################################################################
#
# Initialize drought legacy
#
###############################################################################
"""
    inititialize_legacy!(
                hs::Union{LeafHydraulics,RootHydraulics,StemHydraulics})
    inititialize_legacy!(tree::GrassLikeOrganism)
    inititialize_legacy!(tree::PalmLikeOrganism)
    inititialize_legacy!(tree::TreeLikeOrganism)
    inititialize_legacy!(tree::TreeSimple)

Initialize the drought legacy effects in the xylem, given
- `hs` [`AbstractHydraulicOrgan`] type struct
- `tree` [`TreeSimple`](@ref) type struct
"""
function inititialize_legacy!(
            hs::Union{LeafHydraulics,RootHydraulics,StemHydraulics}
)
    hs.k_history .= 1;
    hs.p_history .= 0;

    return nothing
end




function inititialize_legacy!(tree::GrassLikeOrganism)
    for root in tree.roots
        inititialize_legacy!(root);
    end
    for leaf in tree.leaves
        inititialize_legacy!(leaf);
    end

    return nothing
end




function inititialize_legacy!(tree::PalmLikeOrganism)
    for root in tree.roots
        inititialize_legacy!(root);
    end
    inititialize_legacy!(tree.trunk);
    for leaf in tree.leaves
        inititialize_legacy!(leaf);
    end

    return nothing
end




function inititialize_legacy!(tree::TreeLikeOrganism)
    for root in tree.roots
        inititialize_legacy!(root);
    end
    inititialize_legacy!(tree.trunk);
    for stem in tree.branch
        inititialize_legacy!(stem);
    end
    for leaf in tree.leaves
        inititialize_legacy!(leaf);
    end

    return nothing
end




function inititialize_legacy!(tree::TreeSimple)
    inititialize_legacy!(tree.root);
    inititialize_legacy!(tree.stem);
    inititialize_legacy!(tree.leaf);

    return nothing
end
