#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-May-25: migrate inititialize_legacy! function to clear_legacy!
#
#######################################################################################################################################################################################################
"""
This function clear the legacy for the hydraulic system. Supported methods are

$(METHODLIST)

"""
function clear_legacy! end


#######################################################################################################################################################################################################
#
# Changes to the method
# General
#     2022-May-26: add method for each hydraulic system like leaf, root, and stem
#
#######################################################################################################################################################################################################
"""

    clear_legacy!(hs::Union{LeafHydraulics{FT}, RootHydraulics{FT}, StemHydraulics{FT}}) where {FT<:AbstractFloat}

Clear the legacy for hydraulic system, given
- `hs` `LeafHydraulics`, `RootHydraulics`, or `StemHydraulics` type structure
"""
clear_legacy!(hs::Union{LeafHydraulics{FT}, RootHydraulics{FT}, StemHydraulics{FT}}) where {FT<:AbstractFloat} = (
    hs.k_history .= 1;
    hs.p_history .= 0;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to the method
# General
#     2022-May-26: add method for leaf (hydraulic system nested within)
#
#######################################################################################################################################################################################################
"""

    clear_legacy!(leaf::Leaf{FT}) where {FT<:AbstractFloat}

Clear the legacy for hydraulic system, given
- `leaf` `Leaf` type structure
"""
clear_legacy!(leaf::Leaf{FT}) where {FT<:AbstractFloat} = clear_legacy!(leaf.HS);


#######################################################################################################################################################################################################
#
# Changes to the method
# General
#     2022-May-26: add method for SPAC system (hydraulic system nested within)
#
#######################################################################################################################################################################################################
"""

    clear_legacy!(spac::MonoElementSAPC{FT}) where {FT<:AbstractFloat}

Clear the legacy for SPAC system, given
- `spac` `MonoElementSAPC` type structure
"""
clear_legacy!(spac::MonoElementSAPC{FT}) where {FT<:AbstractFloat} = (
    clear_legacy!(spac.LEAF);
    clear_legacy!(spac.ROOT);
    clear_legacy!(spac.STEM);

    return nothing
)


#######################################################################################################################################################################################################
#
# Changes to the method
# General
#     2022-May-26: add method for SPAC system (hydraulic system nested within)
#
#######################################################################################################################################################################################################
"""

    clear_legacy!(spac::MonoGrassSAPC{FT}) where {FT<:AbstractFloat}

Clear the legacy for SPAC system, given
- `spac` `MonoGrassSAPC` type structure
"""
clear_legacy!(spac::MonoGrassSAPC{FT}) where {FT<:AbstractFloat} = (
    clear_legacy!.(spac.ROOTS);
    clear_legacy!.(spac.LEAVES);

    return nothing
)


#######################################################################################################################################################################################################
#
# Changes to the method
# General
#     2022-May-26: add method for SPAC system (hydraulic system nested within)
#
#######################################################################################################################################################################################################
"""

    clear_legacy!(spac::MonoPalmSAPC{FT}) where {FT<:AbstractFloat}

Clear the legacy for SPAC system, given
- `spac` `MonoPalmSAPC` type structure
"""
clear_legacy!(spac::MonoPalmSAPC{FT}) where {FT<:AbstractFloat} = (
    clear_legacy!.(spac.ROOTS);
    clear_legacy!(spac.TRUNK);
    clear_legacy!.(spac.LEAVES);

    return nothing
)


#######################################################################################################################################################################################################
#
# Changes to the method
# General
#     2022-May-26: add method for SPAC system (hydraulic system nested within)
#
#######################################################################################################################################################################################################
"""

    clear_legacy!(spac::MonoTreeSAPC{FT}) where {FT<:AbstractFloat}

Clear the legacy for SPAC system, given
- `spac` `MonoTreeSAPC` type structure
"""
clear_legacy!(spac::MonoTreeSAPC{FT}) where {FT<:AbstractFloat} = (
    clear_legacy!.(spac.ROOTS);
    clear_legacy!(spac.TRUNK);
    clear_legacy!.(spac.BRANCHES);
    clear_legacy!.(spac.LEAVES);

    return nothing
)
