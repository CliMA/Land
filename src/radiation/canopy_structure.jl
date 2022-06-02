#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jun-02: add abstract type for canopy structure
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractCanopyStructure:
- [`HyperspectralMLCanopy`](@ref)
"""
abstract type AbstractCanopyStructure{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-02: migrate from CanopyLayers
#     2022-Jun-02: rename Canopy4RT to HyperspectralMLCanopy
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save multiple layer hyperspectral canopy parameters

# Fields

$(TYPEDFIELDS)

"""
mutable struct HyperspectralMLCanopy{FT} <: AbstractCanopyStructure{FT}
    # parameters that do not change with time
    "Hot spot parameter"
    HOT_SPOT::FT
    "Leaf inclination angle distribution function parameter a"
    LIDF_A::FT
    "Leaf inclination angle distribution function parameter b"
    LIDF_B::FT
    "Number of azimuth angles"
    N_AZI::Int
    "Number of inclination angles"
    N_INCL::Int
    "Number of canopy layers"
    N_LAYER::Int
    "Inclination angle distribution"
    P_INC::Vector{FT}
    "Clumping structure a"
    Ω_A::FT
    "Clumping structure b"
    Ω_B::FT
    "Mean azimuth angles `[°]`"
    Θ_AZI::Vector{FT}
    "Mean inclination angles `[°]`"
    Θ_INCL::Vector{FT}
    "Bounds of inclination angles `[°]`"
    Θ_INCL_BNDS::Matrix{FT}

    # prognostic variables that change with time
    "Clumping index"
    ci::FT
    "Leaf area index"
    lai::FT

    # caches to speed up calculations
    "Cosine of Θ_AZI"
    _COS_Θ_AZI::Vector{FT}
    "Cosine of Θ_INCL"
    _COS_Θ_INCL::Vector{FT}
    "Sine of Θ_INCL"
    _SIN_Θ_INCL::Vector{FT}
    "Cosine of Θ_AZI - raa"
    _cos_θ_azi_raa::Vector{FT}
    "Cache for volume scatter function"
    _vol_scatter::Vector{FT}
    "Cache for level boundary locations"
    _x_bnds::Vector{FT}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jan-14: add C3 and C4 constructors
#     2022-Jan-24: add C3Cytochrome constructor
#     2022-Jan-24: add p_CO₂_s to the constructor
#     2022-Jan-24: add documentation
#     2022-Feb-07: remove fluorescence model from Leaf struct
#     2022-Feb-11: set default APAR = 1000
#     2022-Feb-11: add colimit option in constructor to enable quick deployment of quadratic colimitation
#     2022-May-25: add leaf hydraulic system into the constructor
#     2022-May-31: add steady state mode option to input options
#     2022-May-25: add new field WIDTH
#
#######################################################################################################################################################################################################
"""

    HyperspectralMLCanopy{FT}(; lai::Number = 3, n_layer::Int = 20, θ_incl_bnds::Matrix = [collect(0:10:80) collect(10:10:90)]) where {FT<:AbstractFloat}

Construct a multiple layer canopy for hyperspectral radiative transfer, given
- `lai` Leaf area index
- `n_layer` Total canopy layers
- `θ_incl_bnds` Inclination angle boundary values
"""
HyperspectralMLCanopy{FT}(; lai::Number = 3, n_layer::Int = 20, θ_incl_bnds::Matrix = [collect(0:10:80) collect(10:10:90)]) where {FT<:AbstractFloat} = (
    _n_incl = size(θ_incl_bnds,1);
    _θ_incl = FT[(θ_incl_bnds[_i,1] + θ_incl_bnds[_i,2]) / 2 for _i in 1:_n_incl];
    _p_incl = ones(_n_incl) / _n_incl;
    _θ_azi  = collect(FT,5:10:360);
    _x_bnds = collect(FT,0:-1/n_layer:-1-eps(FT));

    return HyperspectralMLCanopy{FT}(
                0.05,           # HOT_SPOT
                0,              # LIDF_A
                0,              # LIDF_B
                36,             # N_AZI
                _n_incl,        # N_INCL
                n_layer,        # N_LAYER
                _p_incl,        # P_INC
                1,              # Ω_A
                0,              # Ω_B
                _θ_azi,         # Θ_AZI
                _θ_incl,        # Θ_INCL
                θ_incl_bnds,    # Θ_INCL_BNDS
                1,              # ci
                lai,            # lai
                cosd.(_θ_azi),  # _COS_Θ_AZI
                cosd.(_θ_incl), # _COS_Θ_INCL
                sind.(_θ_incl), # _SIN_Θ_INCL
                cosd.(_θ_azi),  # _cos_θ_azi_raa
                ones(FT,4),     # _vol_scatter
                _x_bnds         # _x_bnds
    )
);
