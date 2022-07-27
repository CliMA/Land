#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-27: add function to initialize SPAC
#     2022-Jun-27: add leaf area controller to make sure soil and leaf areas are consistent with leaf area index
#
#######################################################################################################################################################################################################
"""

    initialize!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat}

Initialize the SPAC, given
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` SPAC

"""
function initialize! end

initialize!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat} = (
    @unpack CANOPY, DIM_LAYER, LEAVES, SOIL = spac;

    # make sure leaf area index setup is correct
    for _i in 1:DIM_LAYER
        LEAVES[_i].HS.AREA = SOIL.AREA * CANOPY.lai / DIM_LAYER;
    end;

    # initialize leaf level spectra
    leaf_spectra!(spac);

    # initialize stomatal conductance
    stomatal_conductance!(spac, FT(0));

    return nothing
);
