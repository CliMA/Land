#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Sep-07: add function to compute GPP for SPAC
#
#######################################################################################################################################################################################################
"""

    GPP(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat}

Return the gross primary productivity per ground area, given
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` SPAC

"""
function GPP end

GPP(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat} = (
    @unpack CANOPY, DIM_LAYER, LEAVES = spac;

    # compute GPP
    _gpp::FT = 0;
    for _i in 1:DIM_LAYER
        _gpp += CANOPY.OPTICS.p_sunlit[_i] * mean(LEAVES[_i].a_gross_sunlit) + (1 - CANOPY.OPTICS.p_sunlit[_i]) * LEAVES[_i].a_gross_shaded;
    end;
    _gpp *= spac.CANOPY.lai / spac.CANOPY.DIM_LAYER;

    return _gpp
);
