#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jul-19: abstractize the leaf
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Abstract type for leaf

Hierarchy of the `AbstractLeaf`
- [`Leaf`](@ref)
- [`Leaves1D`](@ref)
- [`Leaves2D`](@ref)

"""
abstract type AbstractLeaf{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jan-14: refactor the Leaf structure within BIO, PRC, PSM as fields
#     2022-Jan-24: add p_CO₂_s to the structure
#     2022-Jan-24: add FT control to p_CO₂_i
#     2022-Feb-07: moved FLM to PRC
#     2022-May-25: add new field HS, WIDTH
#     2022-Jun-14: use Union instead of Abstract... for type definition
#     2022-Jun-15: add support to BroadbandLeafBiophysics and HyperspectralLeafBiophysics types
#     2022-Jun-29: add APAR_CAR as a field
#     2022-Jun-30: add SM as a field
#     2022-Jul-01: add fields: G_LIMITS, a_gross and a_net
#     2022-Jul-12: add field: ∂g∂t
#     2022-Jul-14: add field: CP, e, cp, and ∂e∂t
#     2022-Jul-19: use kwdef for the constructor
#     2022-Jul-19: remove field p_H₂O_sat
# To do
#     TODO: link leaf water content to BIO_PHYSICS.l_H₂O
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save leaf parameters. This structure is meant for leaf level research and canopy radiative transfer scheme without sunlit and shaded partitioning (ppar and ppar-dependent variables).

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Leaf{FT<:AbstractFloat} <: AbstractLeaf{FT}
    # parameters that do not change with time
    "Whether APAR absorbed by carotenoid is counted as PPAR"
    APAR_CAR::Bool = true
    "[`AbstractLeafBiophysics`](@ref) type leaf biophysical parameters"
    BIO::Union{BroadbandLeafBiophysics{FT}, HyperspectralLeafBiophysics{FT}} = HyperspectralLeafBiophysics{FT}()
    "Specific heat capacity of leaf `[J K⁻¹ kg⁻¹]`"
    CP::FT = 1780
    "Minimal and maximum stomatal conductance for H₂O at 25 °C `[mol m⁻² s⁻¹]`"
    G_LIMITS::Vector{FT} = FT[0.01, 0.3]
    "[`LeafHydraulics`](@ref) type leaf hydraulic system"
    HS::LeafHydraulics{FT} = LeafHydraulics{FT}()
    "[`AbstractReactionCenter`](@ref) type photosynthesis reaction center"
    PRC::Union{VJPReactionCenter{FT}, CytochromeReactionCenter{FT}} = VJPReactionCenter{FT}()
    "[`AbstractPhotosynthesisModel`](@ref) type photosynthesis model"
    PSM::Union{C3VJPModel{FT}, C4VJPModel{FT}, C3CytochromeModel{FT}} = C3VJPModel{FT}()
    "Stomatal model"
    SM::AbstractStomataModel{FT} = WangSM{FT}()
    "Leaf width `[m]`"
    WIDTH::FT = 0.05

    # prognostic variables that change with time
    "Total stored energy per area `[J m⁻²]`"
    e::FT = (CP * BIO.lma * 10 + HS.v_storage * CP_L_MOL(FT)) * T_25()
    "Stomatal conductance to water vapor `[mol m⁻² s⁻¹]`"
    g_H₂O_s::FT = 0.01
    "Absorbed photosynthetically active radiation used for photosynthesis `[μmol m⁻² s⁻¹]`"
    ppar::FT = 1000
    "Current leaf temperature"
    t::FT = T_25()

    # dignostic variables that change with time
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_gross::FT = 0
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_net::FT = 0
    "Combined specific heat capacity of leaf per area `[J K⁻¹ m⁻²]`"
    cp::FT = 0
    "Total leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_CO₂::FT = 0
    "Boundary leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_CO₂_b::FT = 3
    "Leaf internal CO₂ partial pressure `[Pa]`"
    p_CO₂_i::FT = 0
    "Leaf surface CO₂ partial pressure `[Pa]`"
    p_CO₂_s::FT = 0
    "Marginal increase in energy `[W m⁻²]`"
    ∂e∂t::FT = 0
    "Marginal increase of conductance per time `[mol m⁻² s⁻²]`"
    ∂g∂t::FT = 0
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jul-19: clean the constructor function
#
#######################################################################################################################################################################################################
"""

    Leaf{FT}(psm::String, wls::WaveLengthSet{FT} = WaveLengthSet{FT}(); broadband::Bool = false) where {FT<:AbstractFloat}

Constructor for `Leaf`, given
- `psm` Photosynthesis model type, must be `C3`, `C3Cytochrome`, or `C4`
- `wls` [`WaveLengthSet`](@ref) type structure that determines the dimensions of leaf parameters
- `broadband` Whether leaf biophysics is in broadband mode

"""
Leaf{FT}(psm::String, wls::WaveLengthSet{FT} = WaveLengthSet{FT}(); broadband::Bool = false) where {FT<:AbstractFloat} = (
    @assert psm in ["C3", "C3Cytochrome", "C4"] "Photosynthesis model ID must be C3, C4, or C3Cytochrome!";

    _bio = broadband ? HyperspectralLeafBiophysics{FT}(wls) : BroadbandLeafBiophysics{FT}();

    if psm == "C3Cytochrome"
        return Leaf{FT}(BIO = _bio, PRC = CytochromeReactionCenter{FT}(), PSM = C3CytochromeModel{FT}())
    elseif psm == "C4"
        return Leaf{FT}(BIO = _bio, PRC = VJPReactionCenter{FT}(), PSM = C4VJPModel{FT}())
    end;

    return Leaf{FT}(BIO = _bio)
);


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-27: add new structure for leaves with 1D Vector of parameters, such as leaves for sunlit and shaded partitions
#     2022-Jun-27: make BIO BroadbandLeafBiophysics only
#     2022-Jun-28: add a_gross and a_net, make t a Vector, remove _t
#     2022-Jun-30: add a second HS2 for shaded leaves
#     2022-Jun-30: add SM as a field
#     2022-Jul-01: add G_LIMITS as a field
#     2022-Jul-07: make p_H₂O_sat a vector
#     2022-Jul-12: add field: ∂g∂t
#     2022-Jul-14: add field: CP, e, cp, and ∂e∂t
#     2022-Jul-19: use kwdef for the constructor
#     2022-Jul-19: remove field p_H₂O_sat
# To do
#     TODO: link leaf water content to BIO_PHYSICS.l_H₂O
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save leaf parameters for a single canopy layer. This structure is meant for canopy level research and canopy radiative transfer scheme with sunlit and shaded partitioning.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Leaves1D{FT<:AbstractFloat} <: AbstractLeaf{FT}
    # parameters that do not change with time
    "[`BroadbandLeafBiophysics`](@ref) type leaf biophysical parameters"
    BIO::BroadbandLeafBiophysics{FT} = BroadbandLeafBiophysics{FT}()
    "Specific heat capacity of leaf `[J K⁻¹ kg⁻¹]`"
    CP::FT = 1780
    "Minimal and maximum stomatal conductance for H₂O at 25 °C `[mol m⁻² s⁻¹]`"
    G_LIMITS::Vector{FT} = FT[0.01, 0.3]
    "[`LeafHydraulics`](@ref) type leaf hydraulic system"
    HS::LeafHydraulics{FT} = LeafHydraulics{FT}()
    "[`LeafHydraulics`](@ref) type leaf hydraulic system used for other calculations (say sunlit and shaded leaf partitioning)"
    HS2::LeafHydraulics{FT} = LeafHydraulics{FT}()
    "[`AbstractReactionCenter`](@ref) type photosynthesis reaction center"
    PRC::Union{VJPReactionCenter{FT}, CytochromeReactionCenter{FT}} = VJPReactionCenter{FT}()
    "[`AbstractPhotosynthesisModel`](@ref) type photosynthesis model"
    PSM::Union{C3VJPModel{FT}, C4VJPModel{FT}, C3CytochromeModel{FT}} = C3VJPModel{FT}()
    "Stomatal model"
    SM::AbstractStomataModel{FT} = WangSM{FT}()
    "Leaf width `[m]`"
    WIDTH::FT = 0.05

    # prognostic variables that change with time
    "Total stored energy per area `[J m⁻²]`"
    e::Vector{FT} = FT[(CP * BIO.lma * 10 + HS.v_storage * CP_L_MOL(FT)) * T_25(), (CP * BIO.lma * 10 + HS2.v_storage * CP_L_MOL(FT)) * T_25()]
    "Stomatal conductance to water vapor `[mol m⁻² s⁻¹]`"
    g_H₂O_s::Vector{FT} = FT[0.01, 0.01]
    "Absorbed photosynthetically active radiation used for photosynthesis `[μmol m⁻² s⁻¹]`"
    ppar::Vector{FT} = FT[1000, 200]
    "Current leaf temperature"
    t::Vector{FT} = FT[T_25(), T_25()]

    # dignostic variables that change with time
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_gross::Vector{FT} = FT[0, 0]
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_net::Vector{FT} = FT[0, 0]
    "Combined specific heat capacity of leaf per area `[J K⁻¹ m⁻²]`"
    cp::Vector{FT} = FT[0, 0]
    "Total leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_CO₂::Vector{FT} = FT[0, 0]
    "Boundary leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_CO₂_b::Vector{FT} = FT[3, 3]
    "Leaf internal CO₂ partial pressure `[Pa]`"
    p_CO₂_i::Vector{FT} = FT[0, 0]
    "Leaf surface CO₂ partial pressure `[Pa]`"
    p_CO₂_s::Vector{FT} = FT[0, 0]
    "Marginal increase in energy `[W m⁻²]`"
    ∂e∂t::Vector{FT} = FT[0, 0]
    "Marginal increase of conductance per time `[mol m⁻² s⁻²]`"
    ∂g∂t::Vector{FT} = FT[0, 0]
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jul-19: clean the constructor function
#
#######################################################################################################################################################################################################
"""

    Leaves1D{FT}(psm::String) where {FT<:AbstractFloat}

Constructor for `Leaves1D`, given
- `psm` Photosynthesis model type, must be `C3`, `C3Cytochrome`, or `C4`

"""
Leaves1D{FT}(psm::String) where {FT<:AbstractFloat} = (
    @assert psm in ["C3", "C3Cytochrome", "C4"] "Photosynthesis model ID must be C3, C4, or C3Cytochrome!";

    if psm == "C3Cytochrome"
        return Leaves1D{FT}(PRC = CytochromeReactionCenter{FT}(), PSM = C3CytochromeModel{FT}())
    elseif psm == "C4"
        return Leaves1D{FT}(PRC = VJPReactionCenter{FT}(), PSM = C4VJPModel{FT}())
    end;

    return Leaves1D{FT}()
);


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-27: add new structure for leaves with 2D Matrix of parameters for sunlit partitioning and point value for shaded partitioning
#     2022-Jun-27: make BIO HyperspectralLeafBiophysics only
#     2022-Jun-27: add sunlit and shaded ppar to struct (remove the ppar in canopy radiation)
#     2022-Jun-28: add a_gross, a_net, and ϕ_f for sunlit and shaded leaves
#     2022-Jun-29: add APAR_CAR as a field
#     2022-Jun-30: fix documentation
#     2022-Jun-30: add SM as a field
#     2022-Jul-01: add G_LIMITS as a field
#     2022-Jul-12: add fields: ∂g∂t_shaded and ∂g∂t_sunlit
#     2022-Jul-14: add field: CP, e, cp, and ∂e∂t
#     2022-Jul-19: use kwdef for the constructor
#     2022-Jul-19: remove field p_H₂O_sat
# To do
#     TODO: link leaf water content to BIO_PHYSICS.l_H₂O
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save leaf parameters for a single canopy layer. This structure is meant for canopy level research and canopy radiative transfer scheme with sunlit and shaded partitioning as well as leaf
    angular distribution.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Leaves2D{FT<:AbstractFloat} <: AbstractLeaf{FT}
    # parameters that do not change with time
    "Whether APAR absorbed by carotenoid is counted as PPAR"
    APAR_CAR::Bool = true
    "[`HyperspectralLeafBiophysics`](@ref) type leaf biophysical parameters"
    BIO::HyperspectralLeafBiophysics{FT} = HyperspectralLeafBiophysics{FT}()
    "Specific heat capacity of leaf `[J K⁻¹ kg⁻¹]`"
    CP::FT = 1780
    "Minimal and maximum stomatal conductance for H₂O at 25 °C `[mol m⁻² s⁻¹]`"
    G_LIMITS::Vector{FT} = FT[0.01, 0.3]
    "[`LeafHydraulics`](@ref) type leaf hydraulic system"
    HS::LeafHydraulics{FT} = LeafHydraulics{FT}()
    "[`AbstractReactionCenter`](@ref) type photosynthesis reaction center"
    PRC::Union{VJPReactionCenter{FT}, CytochromeReactionCenter{FT}} = VJPReactionCenter{FT}()
    "[`AbstractPhotosynthesisModel`](@ref) type photosynthesis model"
    PSM::Union{C3VJPModel{FT}, C4VJPModel{FT}, C3CytochromeModel{FT}} = C3VJPModel{FT}()
    "Stomatal model"
    SM::AbstractStomataModel{FT} = WangSM{FT}()
    "Leaf width `[m]`"
    WIDTH::FT = 0.05

    # prognostic variables that change with time
    "Total stored energy per area `[J m⁻²]`"
    e::FT = (CP * BIO.lma * 10 + HS.v_storage * CP_L_MOL(FT)) * T_25()
    "Stomatal conductance to water vapor for shaded leaves `[mol m⁻² s⁻¹]`"
    g_H₂O_s_shaded::FT = 0.01
    "Stomatal conductance to water vapor for sunlit leaves `[mol m⁻² s⁻¹]`"
    g_H₂O_s_sunlit::Matrix{FT} = FT(0.01) .* ones(FT, 9, 36)
    "Absorbed photosynthetically active radiation used for photosynthesis for shaded leaves `[μmol m⁻² s⁻¹]`"
    ppar_shaded::FT = 200
    "Absorbed photosynthetically active radiation used for photosynthesis for sunlit leaves `[μmol m⁻² s⁻¹]`"
    ppar_sunlit::Matrix{FT} = 1000 .* ones(FT, 9, 36)
    "Current leaf temperature"
    t::FT = T_25()

    # dignostic variables that change with time
    "Gross photosynthetic rate for shaded leaves `[μmol m⁻² s⁻¹]`"
    a_gross_shaded::FT = 0
    "Gross photosynthetic rate for sunlit leaves `[μmol m⁻² s⁻¹]`"
    a_gross_sunlit::Matrix{FT} = zeros(FT, 9, 36)
    "Net photosynthetic rate for shaded leaves `[μmol m⁻² s⁻¹]`"
    a_net_shaded::FT = 0
    "Net photosynthetic rate for sunlit leaves `[μmol m⁻² s⁻¹]`"
    a_net_sunlit::Matrix{FT} = zeros(FT, 9, 36)
    "Combined specific heat capacity of leaf per area `[J K⁻¹ m⁻²]`"
    cp::FT = 0
    "Total leaf diffusive conductance to CO₂ for shaded leaves `[mol m⁻² s⁻¹]`"
    g_CO₂_shaded::FT = 0
    "Total leaf diffusive conductance to CO₂ for sunlit leaves `[mol m⁻² s⁻¹]`"
    g_CO₂_sunlit::Matrix{FT} = zeros(FT, 9, 36)
    "Boundary leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_CO₂_b::FT = 3
    "Leaf internal CO₂ partial pressure for shaded leaves `[Pa]`"
    p_CO₂_i_shaded::FT = 0
    "Leaf internal CO₂ partial pressure for sunlit leaves `[Pa]`"
    p_CO₂_i_sunlit::Matrix{FT} = zeros(FT, 9, 36)
    "Leaf surface CO₂ partial pressure for shaded leaves `[Pa]`"
    p_CO₂_s_shaded::FT = 0
    "Leaf surface CO₂ partial pressure for sunlit leaves `[Pa]`"
    p_CO₂_s_sunlit::Matrix{FT} = zeros(FT, 9, 36)
    "Fluorescence quantum yield for shaded leaves `[-]`"
    ϕ_f_shaded::FT = 0
    "Fluorescence quantum yield for sunlit leaves `[-]`"
    ϕ_f_sunlit::Matrix{FT} = zeros(FT, 9, 36)
    "Marginal increase in energy `[W m⁻²]`"
    ∂e∂t::FT = 0
    "Marginal increase of conductance per time for shaded leaves `[mol m⁻² s⁻²]`"
    ∂g∂t_shaded::FT = 0
    "Marginal increase of conductance per timefor sunlit leaves `[mol m⁻² s⁻²]`"
    ∂g∂t_sunlit::Matrix{FT} = zeros(FT, 9, 36)

    # caches to speed up calculations
    "Last leaf temperature. If different from t, then make temperature correction"
    _t::FT = 0
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-27: add constructor for Leaves2D
#     2022-Jul-19: clean the constructor function
#
#######################################################################################################################################################################################################
"""

    Leaves2D{FT}(psm::String, wls::WaveLengthSet{FT} = WaveLengthSet{FT}(); n_azi::Int = 36, n_incl::Int = 9) where {FT<:AbstractFloat}

Constructor for `Leaves2D`, given
- `psm` Photosynthesis model type, must be `C3`, `C3Cytochrome`, or `C4`
- `wls` [`WaveLengthSet`](@ref) type structure that determines the dimensions of leaf parameters
- `n_azi` Number of azimuth angles
- `n_incl` Number of inclination angles

"""
Leaves2D{FT}(psm::String, wls::WaveLengthSet{FT} = WaveLengthSet{FT}(); n_azi::Int = 36, n_incl::Int = 9) where {FT<:AbstractFloat} = (
    @assert psm in ["C3", "C3Cytochrome", "C4"] "Photosynthesis model ID must be C3, C4, or C3Cytochrome!";

    if psm == "C3"
        _prc = VJPReactionCenter{FT}();
        _psm = C3VJPModel{FT}();
    elseif psm == "C3Cytochrome"
        _prc = CytochromeReactionCenter{FT}();
        _psm = C3CytochromeModel{FT}();
    elseif psm == "C4"
        _prc = VJPReactionCenter{FT}();
        _psm = C4VJPModel{FT}();
    end;

    return Leaves2D{FT}(
                BIO            = HyperspectralLeafBiophysics{FT}(wls),
                PRC            = _prc,
                PSM            = _psm,
                g_H₂O_s_sunlit = zeros(FT,n_incl,n_azi) .* FT(0.01),
                ppar_sunlit    = zeros(FT,n_incl,n_azi) .* 1000,
                a_gross_sunlit = zeros(FT,n_incl,n_azi),
                a_net_sunlit   = zeros(FT,n_incl,n_azi),
                g_CO₂_sunlit   = zeros(FT,n_incl,n_azi) .* FT(0.01),
                p_CO₂_i_sunlit = zeros(FT,n_incl,n_azi) .+ 20,
                p_CO₂_s_sunlit = zeros(FT,n_incl,n_azi) .+ 40,
                ϕ_f_sunlit     = zeros(FT,n_incl,n_azi),
                ∂g∂t_sunlit    = zeros(FT,n_incl,n_azi)
    )
);
