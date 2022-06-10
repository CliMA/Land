#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-10: migrate the function from CanopyLayers
#     2022-Jun-10: rename function to canopy_fluorescence!
#
#######################################################################################################################################################################################################
"""
This function updates canopy fluorescence profiles. The supported methods include

$(METHODLIST)

"""
function canopy_fluorescence! end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-10: migrate the function from CanopyLayers
#
#######################################################################################################################################################################################################
"""

    canopy_fluorescence!(can::HyperspectralMLCanopy{FT}) where {FT<:AbstractFloat}

Updates canopy radiation profiles for shortwave radiation, given
- `can` `HyperspectralMLCanopy` type struct
"""
canopy_fluorescence!(can::HyperspectralMLCanopy{FT}) where {FT<:AbstractFloat} = (
    ;

    return nothing
);
