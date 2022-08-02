#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: add colimit function back
#     2022-Jan-24: use colimit from psm to abstractize the MinimumColimit and QuadraticColimit methods
#     2022-Feb-07: add C3CytochromeModel support
#     2022-Feb-11: add a wrapper function to step through COLIMIT_CJ and COLIMIT_IP
#     2022-Jul-01: add β to variable list to account for Vmax downregulation used in CLM5
#
#######################################################################################################################################################################################################
"""

    colimit_photosynthesis!(psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}, C4VJPModel{FT}}; β::FT = FT(1)) where {FT<:AbstractFloat}

Colimit the photosynthesis by rubisco-, light-, and product-limited photosynthetic rates, given
- `psm` `C3CytochromeModel`, `C3VJPModel`, or `C4VJPModel` type photosynthesis model
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd (default is 1)

"""
function colimit_photosynthesis! end

colimit_photosynthesis!(psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}, C4VJPModel{FT}}; β::FT = FT(1)) where {FT<:AbstractFloat} =
    colimit_photosynthesis!(psm, psm.COLIMIT_CJ, psm.COLIMIT_IP; β = β);

colimit_photosynthesis!(
            psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}, C4VJPModel{FT}},
            colim_cj::Union{MinimumColimit{FT}, QuadraticColimit{FT}},
            colim_ip::Union{MinimumColimit{FT}, QuadraticColimit{FT}};
            β::FT = FT(1)
) where {FT<:AbstractFloat} = (
    _a_i        = colimited_rate(psm._a_c, psm._a_j, colim_cj);
    psm.a_gross = colimited_rate(psm._a_p, _a_i, colim_ip);
    psm.a_net   = psm.a_gross - β * psm._r_d;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: add colimit function back
#     2022-Feb-07: add C3Cytochrome to method (colimit j_p680 and j_p700 as well)
#     2022-Feb-11: add colimited_rate for general purpose in ETR as well as a_gross
#     2022-Mar-01: add colimit method for serial colimitation
#     2022-Mar-01: add colimit method for square colimitation
#
#######################################################################################################################################################################################################
"""

    colimited_rate(a_1::FT, a_2::FT, colim::MinimumColimit{FT}) where {FT<:AbstractFloat}
    colimited_rate(a_1::FT, a_2::FT, colim::QuadraticColimit{FT}) where {FT<:AbstractFloat}
    colimited_rate(a_1::FT, a_2::FT, colim::SerialColimit{FT}) where {FT<:AbstractFloat}
    colimited_rate(a_1::FT, a_2::FT, colim::SquareColimit{FT}) where {FT<:AbstractFloat}

Return the minimum of two rates, given
- `a_1` Rate 1
- `a_2` Rate 2
- `colim` `MinimumColimit`, `QuadraticColimit`, or `SerialColimit` type struct

"""
function colimited_rate end

colimited_rate(a_1::FT, a_2::FT, colim::MinimumColimit{FT}) where {FT<:AbstractFloat} = min(a_1, a_2);

colimited_rate(a_1::FT, a_2::FT, colim::QuadraticColimit{FT}) where {FT<:AbstractFloat} = lower_quadratic(colim.CURVATURE, -(a_1 + a_2), a_1 * a_2);

colimited_rate(a_1::FT, a_2::FT, colim::SerialColimit{FT}) where {FT<:AbstractFloat} = a_1 * a_2 / (a_1 + a_2);

colimited_rate(a_1::FT, a_2::FT, colim::SquareColimit{FT}) where {FT<:AbstractFloat} = a_1 * a_2 / sqrt(a_1 ^ 2 + a_2 ^ 2);
