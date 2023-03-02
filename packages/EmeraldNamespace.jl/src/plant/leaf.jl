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
#     2022-Feb-07: moved FLM to PRC
#     2022-May-25: add new field HS, WIDTH
#     2022-Jun-14: use Union instead of Abstract... for type definition
#     2022-Jun-15: add support to BroadbandLeafBiophysics and HyperspectralLeafBiophysics types
#     2022-Jun-29: add APAR_CAR as a field
#     2022-Jun-30: add SM as a field
#     2022-Jul-01: add fields: G_LIMITS, a_gross and a_net
#     2022-Jul-12: add field: ∂g∂t
#     2022-Jul-14: add field: CP, e, cp, and ∂e∂t
#     2022-Jul-19: remove field p_H₂O_sat
#     2022-Jul-28: move field _t to PSM
#     2022-Nov-18: use Union type for SM
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save leaf parameters. This structure is meant for leaf level research and canopy radiative transfer scheme without sunlit and shaded partitioning (ppar and ppar-dependent variables).

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Leaf{FT<:AbstractFloat} <: AbstractLeaf{FT}
    # General model information
    "Whether APAR absorbed by carotenoid is counted as PPAR"
    APAR_CAR::Bool = true

    # Constants
    "Specific heat capacity of leaf `[J K⁻¹ kg⁻¹]`"
    CP::FT = 1780
    "Minimal and maximum stomatal conductance for H₂O at 25 °C `[mol m⁻² s⁻¹]`"
    G_LIMITS::Vector{FT} = FT[1e-4, 0.3]
    "Leaf width `[m]`"
    WIDTH::FT = 0.05

    # Embedded structures
    "[`AbstractLeafBiophysics`](@ref) type leaf biophysical parameters"
    BIO::Union{BroadbandLeafBiophysics{FT}, HyperspectralLeafBiophysics{FT}} = HyperspectralLeafBiophysics{FT}()
    "[`LeafHydraulics`](@ref) type leaf hydraulic system"
    HS::LeafHydraulics{FT} = LeafHydraulics{FT}()
    "[`AbstractReactionCenter`](@ref) type photosynthesis reaction center"
    PRC::Union{VJPReactionCenter{FT}, CytochromeReactionCenter{FT}} = VJPReactionCenter{FT}()
    "[`AbstractPhotosynthesisModel`](@ref) type photosynthesis model"
    PSM::Union{C3VJPModel{FT}, C4VJPModel{FT}, C3CytochromeModel{FT}} = C3VJPModel{FT}()
    "Stomatal model"
    SM::Union{AndereggSM{FT}, BallBerrySM{FT}, EllerSM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}, SperrySM{FT}, WangSM{FT}, Wang2SM{FT}} = WangSM{FT}()

    # Prognostic variables (not used for ∂y∂t)
    "Boundary leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_CO₂_b::FT = 3
    "Absorbed photosynthetically active radiation used for photosynthesis `[μmol m⁻² s⁻¹]`"
    ppar::FT = 1000
    "Current leaf temperature"
    t::FT = T₂₅()

    # Prognostic variables (used for ∂y∂t)
    "Total stored energy per area `[J m⁻²]`"
    e::FT = (CP * BIO.lma * 10 + HS.v_storage * CP_L_MOL(FT)) * t
    "Stomatal conductance to water vapor `[mol m⁻² s⁻¹]`"
    g_H₂O_s::FT = 0.01
    "Marginal increase in energy `[W m⁻²]`"
    ∂e∂t::FT = 0
    "Marginal increase of conductance per time `[mol m⁻² s⁻²]`"
    ∂g∂t::FT = 0

    # Diagnostic variables
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_gross::FT = 0
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_net::FT = 0

    # Cache variables
    "Combined specific heat capacity of leaf per area `[J K⁻¹ m⁻²]`"
    _cp::FT = 0
    "Total leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    _g_CO₂::FT = 0
    "Leaf internal CO₂ partial pressure `[Pa]`"
    _p_CO₂_i::FT = 0
    "Leaf surface CO₂ partial pressure `[Pa]`"
    _p_CO₂_s::FT = 0
end


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
#     2022-Jul-19: remove field p_H₂O_sat
#     2022-Nov-18: use Union type for SM
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save leaf parameters for a single canopy layer. This structure is meant for canopy level research and canopy radiative transfer scheme with sunlit and shaded partitioning.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Leaves1D{FT<:AbstractFloat} <: AbstractLeaf{FT}
    # Constants
    "Specific heat capacity of leaf `[J K⁻¹ kg⁻¹]`"
    CP::FT = 1780
    "Minimal and maximum stomatal conductance for H₂O at 25 °C `[mol m⁻² s⁻¹]`"
    G_LIMITS::Vector{FT} = FT[1e-4, 0.3]
    "Leaf width `[m]`"
    WIDTH::FT = 0.05

    # Embedded structures
    "[`BroadbandLeafBiophysics`](@ref) type leaf biophysical parameters"
    BIO::BroadbandLeafBiophysics{FT} = BroadbandLeafBiophysics{FT}()
    "[`LeafHydraulics`](@ref) type leaf hydraulic system"
    HS::LeafHydraulics{FT} = LeafHydraulics{FT}()
    "[`LeafHydraulics`](@ref) type leaf hydraulic system used for other calculations (say sunlit and shaded leaf partitioning)"
    HS2::LeafHydraulics{FT} = LeafHydraulics{FT}()
    "[`AbstractReactionCenter`](@ref) type photosynthesis reaction center"
    PRC::Union{VJPReactionCenter{FT}, CytochromeReactionCenter{FT}} = VJPReactionCenter{FT}()
    "[`AbstractPhotosynthesisModel`](@ref) type photosynthesis model"
    PSM::Union{C3VJPModel{FT}, C4VJPModel{FT}, C3CytochromeModel{FT}} = C3VJPModel{FT}()
    "Stomatal model"
    SM::Union{AndereggSM{FT}, BallBerrySM{FT}, EllerSM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}, SperrySM{FT}, WangSM{FT}, Wang2SM{FT}} = WangSM{FT}()

    # Prognostic variables (not used for ∂y∂t)
    "Boundary leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_CO₂_b::Vector{FT} = FT[3, 3]
    "Absorbed photosynthetically active radiation used for photosynthesis `[μmol m⁻² s⁻¹]`"
    ppar::Vector{FT} = FT[1000, 200]
    "Current leaf temperature"
    t::Vector{FT} = FT[T₂₅(), T₂₅()]

    # Prognostic variables (used for ∂y∂t)
    "Total stored energy per area `[J m⁻²]`"
    e::Vector{FT} = FT[(CP * BIO.lma * 10 + HS.v_storage * CP_L_MOL(FT)) * t[1], (CP * BIO.lma * 10 + HS2.v_storage * CP_L_MOL(FT)) * t[2]]
    "Stomatal conductance to water vapor `[mol m⁻² s⁻¹]`"
    g_H₂O_s::Vector{FT} = FT[0.01, 0.01]
    "Marginal increase in energy `[W m⁻²]`"
    ∂e∂t::Vector{FT} = FT[0, 0]
    "Marginal increase of conductance per time `[mol m⁻² s⁻²]`"
    ∂g∂t::Vector{FT} = FT[0, 0]

    # Diagnostic variables
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_gross::Vector{FT} = FT[0, 0]
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_net::Vector{FT} = FT[0, 0]

    # Cache variables
    "Combined specific heat capacity of leaf per area `[J K⁻¹ m⁻²]`"
    _cp::Vector{FT} = FT[0, 0]
    "Total leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    _g_CO₂::Vector{FT} = FT[0, 0]
    "Leaf internal CO₂ partial pressure `[Pa]`"
    _p_CO₂_i::Vector{FT} = FT[0, 0]
    "Leaf surface CO₂ partial pressure `[Pa]`"
    _p_CO₂_s::Vector{FT} = FT[0, 0]
end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-27: add new structure for leaves with 2D Matrix of parameters for sunlit partitioning and point value for shaded partitioning
#     2022-Jun-27: make BIO HyperspectralLeafBiophysics only
#     2022-Jun-27: add sunlit and shaded ppar to struct (remove the ppar in canopy radiation)
#     2022-Jun-28: add a_gross, a_net, and ϕ_f for sunlit and shaded leaves
#     2022-Jun-29: add APAR_CAR as a field
#     2022-Jun-30: add SM as a field
#     2022-Jul-01: add G_LIMITS as a field
#     2022-Jul-12: add fields: ∂g∂t_shaded and ∂g∂t_sunlit
#     2022-Jul-14: add field: CP, e, cp, and ∂e∂t
#     2022-Jul-19: remove field p_H₂O_sat
#     2022-Jul-19: add dimension control to struct
#     2022-Jul-28: move field _t to PSM
#     2022-Nov-18: use Union type for SM
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
    # Dimensions
    "Dimension of azimuth angles"
    DIM_AZI::Int = 36
    "Dimension of inclination angles"
    DIM_INCL::Int = 9

    # General model information
    "Whether APAR absorbed by carotenoid is counted as PPAR"
    APAR_CAR::Bool = true

    # Constants
    "Specific heat capacity of leaf `[J K⁻¹ kg⁻¹]`"
    CP::FT = 1780
    "Minimal and maximum stomatal conductance for H₂O at 25 °C `[mol m⁻² s⁻¹]`"
    G_LIMITS::Vector{FT} = FT[1e-4, 0.3]
    "Leaf width `[m]`"
    WIDTH::FT = 0.05

    # Embedded structures
    "[`HyperspectralLeafBiophysics`](@ref) type leaf biophysical parameters"
    BIO::HyperspectralLeafBiophysics{FT} = HyperspectralLeafBiophysics{FT}()
    "[`LeafHydraulics`](@ref) type leaf hydraulic system"
    HS::LeafHydraulics{FT} = LeafHydraulics{FT}()
    "[`AbstractReactionCenter`](@ref) type photosynthesis reaction center"
    PRC::Union{VJPReactionCenter{FT}, CytochromeReactionCenter{FT}} = VJPReactionCenter{FT}()
    "[`AbstractPhotosynthesisModel`](@ref) type photosynthesis model"
    PSM::Union{C3VJPModel{FT}, C4VJPModel{FT}, C3CytochromeModel{FT}} = C3VJPModel{FT}()
    "Stomatal model"
    SM::Union{AndereggSM{FT}, BallBerrySM{FT}, EllerSM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}, SperrySM{FT}, WangSM{FT}, Wang2SM{FT}} = WangSM{FT}()

    # Prognostic variables (not used for ∂y∂t)
    "Boundary leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_CO₂_b::FT = 3
    "Absorbed photosynthetically active radiation used for photosynthesis for shaded leaves `[μmol m⁻² s⁻¹]`"
    ppar_shaded::FT = 200
    "Absorbed photosynthetically active radiation used for photosynthesis for sunlit leaves `[μmol m⁻² s⁻¹]`"
    ppar_sunlit::Matrix{FT} = 1000 .* ones(FT, DIM_INCL, DIM_AZI)
    "Current leaf temperature `[K]`"
    t::FT = T₂₅()

    # Prognostic variables (used for ∂y∂t)
    "Total stored energy per area `[J m⁻²]`"
    e::FT = (CP * BIO.lma * 10 + HS.v_storage * CP_L_MOL(FT)) * t
    "Stomatal conductance to water vapor for shaded leaves `[mol m⁻² s⁻¹]`"
    g_H₂O_s_shaded::FT = 0.01
    "Stomatal conductance to water vapor for sunlit leaves `[mol m⁻² s⁻¹]`"
    g_H₂O_s_sunlit::Matrix{FT} = FT(0.01) .* ones(FT, DIM_INCL, DIM_AZI)
    "Marginal increase in energy `[W m⁻²]`"
    ∂e∂t::FT = 0
    "Marginal increase of conductance per time for shaded leaves `[mol m⁻² s⁻²]`"
    ∂g∂t_shaded::FT = 0
    "Marginal increase of conductance per timefor sunlit leaves `[mol m⁻² s⁻²]`"
    ∂g∂t_sunlit::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)

    # Diagnostic variables
    "Gross photosynthetic rate for shaded leaves `[μmol m⁻² s⁻¹]`"
    a_gross_shaded::FT = 0
    "Gross photosynthetic rate for sunlit leaves `[μmol m⁻² s⁻¹]`"
    a_gross_sunlit::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
    "Net photosynthetic rate for shaded leaves `[μmol m⁻² s⁻¹]`"
    a_net_shaded::FT = 0
    "Net photosynthetic rate for sunlit leaves `[μmol m⁻² s⁻¹]`"
    a_net_sunlit::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
    "Fluorescence quantum yield for shaded leaves `[-]`"
    ϕ_f_shaded::FT = 0
    "Fluorescence quantum yield for sunlit leaves `[-]`"
    ϕ_f_sunlit::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)

    # Cache variables
    "Combined specific heat capacity of leaf per area `[J K⁻¹ m⁻²]`"
    _cp::FT = 0
    "Total leaf diffusive conductance to CO₂ for shaded leaves `[mol m⁻² s⁻¹]`"
    _g_CO₂_shaded::FT = 0
    "Total leaf diffusive conductance to CO₂ for sunlit leaves `[mol m⁻² s⁻¹]`"
    _g_CO₂_sunlit::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
    "Leaf internal CO₂ partial pressure for shaded leaves `[Pa]`"
    _p_CO₂_i_shaded::FT = 0
    "Leaf internal CO₂ partial pressure for sunlit leaves `[Pa]`"
    _p_CO₂_i_sunlit::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
    "Leaf surface CO₂ partial pressure for shaded leaves `[Pa]`"
    _p_CO₂_s_shaded::FT = 0
    "Leaf surface CO₂ partial pressure for sunlit leaves `[Pa]`"
    _p_CO₂_s_sunlit::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
end
