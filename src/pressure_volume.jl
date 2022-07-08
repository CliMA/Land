#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-May-24: migrate function to version v0.3
#     2022-May-24: rename the function to xylem_pressure
#     2022-May-24: add method for LinearPVCurve
#     2022-May-24: add method for SegmentedPVCurve
#     2022-Jul-08: deflate documentations
#
#######################################################################################################################################################################################################
"""

    xylem_pressure(pv::LinearPVCurve{FT}, rvol::FT, T::FT) where {FT<:AbstractFloat}
    xylem_pressure(pv::SegmentedPVCurve{FT}, rvol::FT, T::FT) where {FT<:AbstractFloat}

Return the xylem water pressure in MPa, given
- `pv` `LinearPVCurve` or `SegmentedPVCurve` type pressure volume curve
- `rvol` Relative volume (relative to maximum)
- `T` Temperature

"""
function xylem_pressure end

xylem_pressure(pv::LinearPVCurve{FT}, rvol::FT, T::FT) where {FT<:AbstractFloat} = (rvol - 1) / pv.SLOPE;

xylem_pressure(pv::SegmentedPVCurve{FT}, rvol::FT, T::FT) where {FT<:AbstractFloat} = (
    @unpack C_ALL, RWC_APO, RWC_TLP, Ε_BULK = pv;

    if rvol > RWC_TLP
        return -C_ALL * GAS_R(FT) * T / (rvol - RWC_APO) * FT(1e-6) + Ε_BULK * (rvol - RWC_TLP)
    elseif rvol > RWC_APO
        return -C_ALL * GAS_R(FT) * T / (rvol - RWC_APO) * FT(1e-6)
    else
        return FT(-100)
    end
);


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-May-27: migrate function to version v0.3
#     2022-May-27: rename the function to capacitance_buffer
#     2022-May-27: add method for LinearPVCurve
#     2022-May-27: add method for SegmentedPVCurve
#     2022-May-31: add documentation
#     2022-Jul-08: deflate documentations
#
#######################################################################################################################################################################################################
"""

    capacitance_buffer(pvc::LinearPVCurve{FT}) where {FT<:AbstractFloat}
    capacitance_buffer(pvc::SegmentedPVCurve{FT}) where {FT<:AbstractFloat}

Return the relative capacictance buffer rate, given
- `pv` `LinearPVCurve` or `SegmentedPVCurve` type pressure volume curve

"""
function capacitance_buffer end

capacitance_buffer(pvc::LinearPVCurve{FT}) where {FT<:AbstractFloat} = pvc.K_REFILL;

capacitance_buffer(pvc::SegmentedPVCurve{FT}) where {FT<:AbstractFloat} = pvc.K_REFILL * (1 - pvc.RWC_APO);
