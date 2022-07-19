#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jun-02: add abstract type for LIDF algorithms
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractLIDFAlgorithm:
- [`VerhoefLIDF`](@ref)
"""
abstract type AbstractLIDFAlgorithm{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-02: migrate from CanopyLayers
#     2022-Jun-02: rename Canopy4RT to HyperspectralMLCanopy
#     2022-Jun-02: abstractize LIDF as a field
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for Verhoef LIDF algorithm

# Fields

$(TYPEDFIELDS)

"""
mutable struct VerhoefLIDF{FT<:AbstractFloat} <: AbstractLIDFAlgorithm{FT}
    # parameters that do not change with time
    "Leaf inclination angle distribution function parameter a"
    A::FT
    "Leaf inclination angle distribution function parameter b"
    B::FT
end


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jun-02: add abstract type for canopy structure
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractCanopy:
- [`BroadbandSLCanopy`](@ref)
- [`HyperspectralMLCanopy`](@ref)
"""
abstract type AbstractCanopy{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-15: add struct for broadband radiative transfer scheme such as two leaf model
#     2022-Jun-15: add more cache variables
#     2022-Jun-15: add radiation profile
#     2022-Jun-15: remove RATIO_HV to compute the coefficient numerically
#     2022-Jun-16: remove some cache variables
#     2022-Jun-16: add fields: Θ_INCL_BNDS
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save single layer broadband canopy parameters

# Fields

$(TYPEDFIELDS)

"""
mutable struct BroadbandSLCanopy{FT<:AbstractFloat} <: AbstractCanopy{FT}
    # parameters that do not change with time
    "Leaf inclination angle distribution function algorithm"
    LIDF::Union{VerhoefLIDF{FT}}
    "Inclination angle distribution"
    P_INCL::Vector{FT}
    "Canopy radiation profiles"
    RADIATION::BroadbandSLCanopyRadiationProfile{FT}
    "Mean inclination angles `[°]`"
    Θ_INCL::Vector{FT}
    "Bounds of inclination angles `[°]`"
    Θ_INCL_BNDS::Matrix{FT}

    # prognostic variables that change with time
    "Clumping index"
    ci::FT
    "Leaf area index"
    lai::FT
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-15: add constructor
#     2022-Jun-15: add more cache variables
#     2022-Jun-15: add radiation profile
#     2022-Jun-15: remove RATIO_HV to compute the coefficient numerically
#     2022-Jun-16: remove some cache variables
#     2022-Jun-16: add fields: Θ_INCL_BNDS
#
#######################################################################################################################################################################################################
"""

    BroadbandSLCanopy{FT}(; lai::Number = 3, θ_incl_bnds::Matrix = [collect(0:10:80) collect(10:10:90)]) where {FT<:AbstractFloat}

Construct a single layer canopy for hyperspectral radiative transfer, given
- `lai` Leaf area index
- `θ_incl_bnds` Inclination angle boundary values
"""
BroadbandSLCanopy{FT}(; lai::Number = 3, θ_incl_bnds::Matrix = [collect(0:10:80) collect(10:10:90)]) where {FT<:AbstractFloat} = (
    _n_incl = size(θ_incl_bnds,1);
    _θ_incl = FT[(θ_incl_bnds[_i,1] + θ_incl_bnds[_i,2]) / 2 for _i in 1:_n_incl];
    _p_incl = ones(_n_incl) / _n_incl;
    _rad = BroadbandSLCanopyRadiationProfile{FT}(n_incl = _n_incl);

    return BroadbandSLCanopy{FT}(
                VerhoefLIDF{FT}(0,0),   # LIDF
                _p_incl,                # P_INCL
                _rad,                   # RADIATION
                _θ_incl,                # Θ_INCL
                θ_incl_bnds,            # Θ_INCL_BNDS
                1,                      # ci
                lai                     # lai
    )
);


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-02: migrate from CanopyLayers
#     2022-Jun-02: rename Canopy4RT to HyperspectralMLCanopy
#     2022-Jun-02: abstractize LIDF as a field
#     2022-Jun-07: add cache variable _1_AZI, _COS²_Θ_INCL, _COS_Θ_INCL_AZI, _COS²_Θ_INCL_AZI
#     2022-Jun-07: remove cache variable _cos_θ_azi_raa, _vol_scatter
#     2022-Jun-09: add new field: APAR_CAR, RADIATION, WLSET
#     2022-Jun-13: use Union instead of Abstract... for type definition
#     2022-Jun-15: rename to HyperspectralMLCanopyOpticalProperty and HyperspectralMLCanopyRadiationProfile
#     2022-Jun-16: remove some cache variables
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save multiple layer hyperspectral canopy parameters

# Fields

$(TYPEDFIELDS)

"""
mutable struct HyperspectralMLCanopy{FT<:AbstractFloat} <: AbstractCanopy{FT}
    # parameters that do not change with time
    "Whether Carotenoid absorption is accounted for in APAR"
    APAR_CAR::Bool
    "Hot spot parameter"
    HOT_SPOT::FT
    "Leaf inclination angle distribution function algorithm"
    LIDF::Union{VerhoefLIDF{FT}}
    "Number of azimuth angles"
    N_AZI::Int
    "Number of inclination angles"
    N_INCL::Int
    "Number of canopy layers"
    N_LAYER::Int
    "Canopy optical properties"
    OPTICS::HyperspectralMLCanopyOpticalProperty{FT}
    "Inclination angle distribution"
    P_INCL::Vector{FT}
    "Canopy radiation profiles"
    RADIATION::HyperspectralMLCanopyRadiationProfile{FT}
    "Wave length set used to paramertize other variables"
    WLSET::WaveLengthSet{FT}
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
    "Ones with the length of Θ_AZI"
    _1_AZI::Vector{FT}
    "Cosine of Θ_AZI"
    _COS_Θ_AZI::Vector{FT}
    "Square of cosine of Θ_INCL"
    _COS²_Θ_INCL::Vector{FT}
    "Square of cosine of Θ_INCL at different azimuth angles"
    _COS²_Θ_INCL_AZI::Matrix{FT}
    "Cache for level boundary locations"
    _x_bnds::Vector{FT}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-02: add constructor
#     2022-Jun-02: abstractize LIDF as a field
#     2022-Jun-07: add cache variable _1_AZI, _COS²_Θ_INCL, _COS_Θ_INCL_AZI, _COS²_Θ_INCL_AZI
#     2022-Jun-07: remove cache variable _cos_θ_azi_raa, _vol_scatter
#     2022-Jun-08: add n_λ to options to initialize CanopyOpticalProperty field
#     2022-Jun-09: add new field: APAR_CAR, RADIATION, WLSET
#     2022-Jun-10: remove n_λ from options and use the N in wls
#     2022-Jun-10: add SIF excitation and fluorescence length control
#     2022-Jun-15: fix documentation
#     2022-Jun-15: rename to HyperspectralMLCanopyOpticalProperty and HyperspectralMLCanopyRadiationProfile
#     2022-Jun-16: remove some cache variables
#
#######################################################################################################################################################################################################
"""

    HyperspectralMLCanopy{FT}(
                wls::WaveLengthSet{FT} = WaveLengthSet{FT}();
                lai::Number = 3,
                n_layer::Int = 20,
                θ_incl_bnds::Matrix = [collect(0:10:80) collect(10:10:90)]
    ) where {FT<:AbstractFloat}

Construct a multiple layer canopy for hyperspectral radiative transfer, given
- `wls` [`WaveLengthSet`](@ref) type struct that defines wavelength settings
- `lai` Leaf area index
- `n_layer` Total canopy layers
- `θ_incl_bnds` Inclination angle boundary values
"""
HyperspectralMLCanopy{FT}(
            wls::WaveLengthSet{FT} = WaveLengthSet{FT}();
            lai::Number = 3,
            n_layer::Int = 20,
            θ_incl_bnds::Matrix = [collect(0:10:80) collect(10:10:90)]
) where {FT<:AbstractFloat} = (
    _n_incl  = size(θ_incl_bnds,1);
    _θ_incl  = FT[(θ_incl_bnds[_i,1] + θ_incl_bnds[_i,2]) / 2 for _i in 1:_n_incl];
    _p_incl  = ones(_n_incl) / _n_incl;
    _θ_azi   = collect(FT,5:10:360);
    _x_bnds  = collect(FT,0:-1/n_layer:-1-eps(FT));
    _can_opt = HyperspectralMLCanopyOpticalProperty{FT}(DIM_AZI = 36, DIM_INCL = _n_incl, DIM_LAYER = n_layer, DIM_SIF = wls.NΛ_SIF, DIM_SIFE = wls.NΛ_SIFE, DIM_WL = wls.NΛ);
    _can_rad = HyperspectralMLCanopyRadiationProfile{FT}(; n_layer = n_layer, n_par = wls.NΛ_PAR, n_λ = wls.NΛ, n_λf = wls.NΛ_SIF);
    _cos_θ   = cosd.(_θ_incl);
    _cos²_θ  = _cos_θ .^ 2;

    return HyperspectralMLCanopy{FT}(
                true,                   # APAR_CAR
                0.05,                   # HOT_SPOT
                VerhoefLIDF{FT}(0,0),   # LIDF
                36,                     # N_AZI
                _n_incl,                # N_INCL
                n_layer,                # N_LAYER
                _can_opt,               # OPTICS
                _p_incl,                # P_INCL
                _can_rad,               # RADIATION
                wls,                    # WLSET
                1,                      # Ω_A
                0,                      # Ω_B
                _θ_azi,                 # Θ_AZI
                _θ_incl,                # Θ_INCL
                θ_incl_bnds,            # Θ_INCL_BNDS
                1,                      # ci
                lai,                    # lai
                ones(FT,36),            # _1_AZI
                cosd.(_θ_azi),          # _COS_Θ_AZI
                _cos²_θ,                # _COS²_Θ_INCL
                _cos²_θ * ones(FT,36)', # _COS²_Θ_INCL_AZI
                _x_bnds                 # _x_bnds
    )
);
