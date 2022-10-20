#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Oct-19: add function to compute canopy net primary productivity
#
#######################################################################################################################################################################################################
"""

    CNPP(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat}

Return the canopy net primary productivity per ground area, given
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` SPAC

"""
function CNPP end

CNPP(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat} = (
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


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Oct-19: add function to compute canopy integrated PPAR
#
#######################################################################################################################################################################################################
"""

    PPAR(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat}

Return the canopy integrated PPAR per ground area, given
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` SPAC

"""
function PPAR end

PPAR(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat} = (
    @unpack CANOPY, DIM_LAYER, LEAVES = spac;

    # compute GPP
    _ppar::FT = 0;
    for _i in 1:DIM_LAYER
        _ppar += CANOPY.OPTICS.p_sunlit[_i] * mean(LEAVES[_i].ppar_sunlit) + (1 - CANOPY.OPTICS.p_sunlit[_i]) * LEAVES[_i].ppar_shaded;
    end;
    _ppar *= spac.CANOPY.lai / spac.CANOPY.DIM_LAYER;

    return _ppar
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Sep-08: add function to add up transpiration rate for SPAC
#
#######################################################################################################################################################################################################
"""

    T_VEG(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat}

Return the transpiration rate per ground area, given
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` SPAC

"""
function T_VEG end

T_VEG(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat} = (
    @unpack CANOPY, DIM_LAYER, LEAVES = spac;

    # compute transpiration rate
    _tran::FT = 0;
    for _i in 1:DIM_LAYER
        _tran += flow_out(LEAVES[_i]);
    end;
    _tran *= spac.CANOPY.lai / spac.CANOPY.DIM_LAYER;

    return _tran
);
