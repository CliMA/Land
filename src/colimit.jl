#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: add colimit function back
#
#######################################################################################################################################################################################################
"""
This function update gross and net photosynthetic rates using different colimitation method. Supported methods are

$(METHODLIST)

"""
function colimit_photosynthesis! end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jan-24: use colimit from psm to abstractize the MinimumColimit and QuadraticColimit methods
#     2022-Jan-24: fix documentation
#     2022-Feb-07: add C3CytochromeModel support
#
#######################################################################################################################################################################################################
"""

    colimit_photosynthesis!(psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}, C4VJPModel{FT}}) where {FT<:AbstractFloat}

Colimit the photosynthesis by rubisco-, light-, and product-limited photosynthetic rates, given
- `psm` `C3VJPModel` or `C4VJPModel` type photosynthesis model
"""
colimit_photosynthesis!(psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}, C4VJPModel{FT}}) where {FT<:AbstractFloat} = (
    colimit_photosynthesis!(psm, psm.COLIMIT_CJ, psm.COLIMIT_IP);

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Feb-11: add a wrapper function to step through COLIMIT_CJ and COLIMIT_IP
#
#######################################################################################################################################################################################################
"""

    colimit_photosynthesis!(
                psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}, C4VJPModel{FT}},
                colim_cj::Union{MinimumColimit{FT}, QuadraticColimit{FT}},
                colim_ip::Union{MinimumColimit{FT}, QuadraticColimit{FT}}
    ) where {FT<:AbstractFloat}

Colimit the photosynthesis by rubisco-, light-, and product-limited photosynthetic rates, given
- `psm` `C3VJPModel` or `C4VJPModel` type photosynthesis model
- `colim_cj` Colimit method for rubisco- and light-limited rates
- `colim_ip` Colimit method for the colimited rubisco- and light-limited rate and product-limited rate
"""
colimit_photosynthesis!(
            psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}, C4VJPModel{FT}},
            colim_cj::Union{MinimumColimit{FT}, QuadraticColimit{FT}},
            colim_ip::Union{MinimumColimit{FT}, QuadraticColimit{FT}}
) where {FT<:AbstractFloat} = (
    _a_i        = colimited_rate(psm.a_c, psm.a_j, colim_cj);
    psm.a_gross = colimited_rate(psm.a_p, _a_i, colim_ip);
    psm.a_net   = psm.a_gross - psm.r_d;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Feb-11: add colimited_rate for general purpose in ETR as well as a_gross
#
#######################################################################################################################################################################################################
"""
This function returns the colimited rates from two rates. Supported methods are

$(METHODLIST)

"""
function colimited_rate end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jan-14: add colimit function back
#     2022-Jan-24: fix documentation
#     2022-Feb-07: add C3Cytochrome to method (colimit j_p680 and j_p700 as well)
#     2022-Feb-11: refactor this function to return the minimum of two rates
#
#######################################################################################################################################################################################################
"""

    colimited_rate(a_1::FT, a_2::FT, colim::MinimumColimit{FT}) where {FT<:AbstractFloat}

Return the minimum of two rates, given
- `a_1` Rate 1
- `a_2` Rate 2
- `colim` `MinimumColimit` type struct
"""
colimited_rate(a_1::FT, a_2::FT, colim::MinimumColimit{FT}) where {FT<:AbstractFloat} = min(a_1, a_2);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jan-14: add colimit function back
#     2022-Jan-14: unpack CONSTANT from the input variables only
#     2022-Jan-24: fix documentation
#     2022-Feb-07: add C3Cytochrome to method (colimit j_p680 and j_p700 as well)
#     2022-Feb-11: refactor this function to return the quadratic colimitation of two rates
#
#######################################################################################################################################################################################################
"""

    colimited_rate(a_1::FT, a_2::FT, colim::QuadraticColimit{FT}) where {FT<:AbstractFloat}

Return the quadratic colimitation of two rates, given
- `a_1` Rate 1
- `a_2` Rate 2
- `colim` `QuadraticColimit` type struct
"""
colimited_rate(a_1::FT, a_2::FT, colim::QuadraticColimit{FT}) where {FT<:AbstractFloat} = lower_quadratic(colim.CURVATURE, -(a_1 + a_2), a_1 * a_2);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Mar-01: add colimit method for serial colimitation
#
#######################################################################################################################################################################################################
"""

    colimited_rate(a_1::FT, a_2::FT, colim::SerialColimit{FT}) where {FT<:AbstractFloat}

Return the serial colimitation of two rates, given
- `a_1` Rate 1
- `a_2` Rate 2
- `colim` `SerialColimit` type struct
"""
colimited_rate(a_1::FT, a_2::FT, colim::SerialColimit{FT}) where {FT<:AbstractFloat} = a_1 * a_2 / (a_1 + a_2);
