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


clear_legacy!(hs::LeafHydraulics{FT}) where {FT<:AbstractFloat} = (
    hs.k_history .= 1;
    hs.p_history .= 0;

    return nothing
);


clear_legacy!(leaf::Leaf{FT}) where {FT<:AbstractFloat} = clear_legacy!(leaf.HS);
