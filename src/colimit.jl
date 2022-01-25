#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: add colimt function back
#
#######################################################################################################################################################################################################
"""
This function update gross and net photosynthetic rates using different colimitation method. Supported methods are

$(METHODLIST)

"""
function colimit_photosynthesis! end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-24: use colimit from psm to abstractize the MinimumColimit and QuadraticColimit methods
#     2022-Jan-24: fix documentation
#
#######################################################################################################################################################################################################
"""

    colimit_photosynthesis!(psm::Union{C3VJPModel{FT}, C4VJPModel{FT}}) where {FT<:AbstractFloat}

Use the minimal photosynthetic rates of rubisco-, light-, and product-limited photosynthetic rates, given
- `psm` `C3VJPModel` or `C4VJPModel` type photosynthesis model
"""
colimit_photosynthesis!(psm::Union{C3VJPModel{FT}, C4VJPModel{FT}}) where {FT<:AbstractFloat} = (
    colimit_photosynthesis!(psm, psm.COLIMIT);

    return nothing;
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: add colimt function back
#     2022-Jan-24: fix documentation
#
#######################################################################################################################################################################################################
"""

    colimit_photosynthesis!(psm::Union{C3VJPModel{FT}, C4VJPModel{FT}}, colim::MinimumColimit{FT}) where {FT<:AbstractFloat}

Use the minimal photosynthetic rates of rubisco-, light-, and product-limited photosynthetic rates, given
- `psm` `C3VJPModel` or `C4VJPModel` type photosynthesis model
- `colim` `MinimumColimit` type struct
"""
colimit_photosynthesis!(psm::Union{C3VJPModel{FT}, C4VJPModel{FT}}, colim::MinimumColimit{FT}) where {FT<:AbstractFloat} = (
    psm.a_gross = min(psm.a_c, psm.a_j, psm.a_p);
    psm.a_net   = psm.a_gross - psm.r_d;

    return nothing;
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: add colimt function back
#     2022-Jan-14: unpack CONSTANT from the input variables only
#     2022-Jan-24: fix documentation
#
#######################################################################################################################################################################################################
"""

    colimit_photosynthesis!(psm::Union{C3VJPModel{FT}, C4VJPModel{FT}}, colim::QuadraticColimit{FT}) where {FT<:AbstractFloat}

Use the minimal photosynthetic rates of rubisco-, light-, and product-limited photosynthetic rates, given
- `psm` `C3VJPModel` or `C4VJPModel` type photosynthesis model
- `colim` `QuadraticColimit` type struct
"""
colimit_photosynthesis!(psm::Union{C3VJPModel{FT}, C4VJPModel{FT}}, colim::QuadraticColimit{FT}) where {FT<:AbstractFloat} = (
    @unpack CURVATURE = colim;

    _a = lower_quadratic(CURVATURE, -(psm.a_c + psm.a_j), psm.a_c * psm.a_j);
    _a = lower_quadratic(CURVATURE, -(_a + psm.a_p), _a  * psm.a_p);

    psm.a_gross = isnan(_a) ? min(psm.a_c, psm.a_j, psm.a_p) : _a;
    psm.a_net   = psm.a_gross - psm.r_d;

    return nothing
);
