#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Oct-19: add function to compute canopy net primary productivity
#
#######################################################################################################################################################################################################
"""

    canopy_net_primary_productivity(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat}

Return the canopy net primary productivity per ground area, given
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` SPAC

"""
function canopy_net_primary_productivity end

canopy_net_primary_productivity(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat} = (
    @unpack CANOPY, DIM_LAYER, LEAVES = spac;

    # compute GPP
    _cnpp::FT = 0;
    for _i in 1:DIM_LAYER
        _cnpp += CANOPY.OPTICS.p_sunlit[_i] * mean(LEAVES[_i].a_net_sunlit) + (1 - CANOPY.OPTICS.p_sunlit[_i]) * LEAVES[_i].a_net_shaded;
    end;
    _cnpp *= spac.CANOPY.lai / spac.CANOPY.DIM_LAYER;

    return _cnpp
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Sep-07: add function to add up GPP for SPAC
#
#######################################################################################################################################################################################################
"""

    gross_primary_productivity(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat}

Return the gross primary productivity per ground area, given
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` SPAC

"""
function gross_primary_productivity end

gross_primary_productivity(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat} = (
    @unpack CANOPY, DIM_LAYER, LEAVES = spac;

    # compute GPP
    _gpp::FT = 0;
    for _i in 1:DIM_LAYER
        _gpp += CANOPY.OPTICS.p_sunlit[_i] * mean(LEAVES[_i].a_gross_sunlit) + (1 - CANOPY.OPTICS.p_sunlit[_i]) * LEAVES[_i].a_gross_shaded;
    end;
    _gpp *= spac.CANOPY.lai / spac.CANOPY.DIM_LAYER;

    return _gpp
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Sep-08: add function to add up transpiration rate for SPAC
#
#######################################################################################################################################################################################################
"""

    transpiration_rate(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat}

Return the transpiration rate per ground area, given
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` SPAC

"""
function transpiration_rate end

transpiration_rate(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat} = (
    @unpack CANOPY, DIM_LAYER, LEAVES = spac;

    # compute transpiration rate
    _tran::FT = 0;
    for _i in 1:DIM_LAYER
        _tran += flow_out(LEAVES[_i]);
    end;
    _tran *= spac.CANOPY.lai / spac.CANOPY.DIM_LAYER;

    return _tran
);
