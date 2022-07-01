#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jan-14: refactor the Leaf structure within BIO, PRC, PSM as fields
#     2022-Jan-24: add p_CO₂_s to the structure
#     2022-Jan-24: fix documentation
#     2022-Feb-07: moved FLM to PRC
#     2022-May-25: add new field HS
#     2022-May-25: add new field WIDTH
#     2022-Jun-14: use Union instead of Abstract... for type definition
#     2022-Jun-15: add support to BroadbandLeafBiophysics and HyperspectralLeafBiophysics types
#     2022-Jun-29: add APAR_CAR as a field
#     2022-Jun-30: add SM as a field
#     2022-Jul-01: add G_LIMITS as a field
# Bug fixes:
#     2022-Jan-24: add FT control to p_CO₂_i
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
mutable struct Leaf{FT<:AbstractFloat}
    # parameters that do not change with time
    "Whether APAR absorbed by carotenoid is counted as PPAR"
    APAR_CAR::Bool
    "[`AbstractLeafBiophysics`](@ref) type leaf biophysical parameters"
    BIO::Union{BroadbandLeafBiophysics{FT}, HyperspectralLeafBiophysics{FT}}
    "Minimal and maximum stomatal conductance for H₂O at 25 °C `[mol m⁻² s⁻¹]`"
    G_LIMITS::Vector{FT}
    "[`LeafHydraulics`](@ref) type leaf hydraulic system"
    HS::LeafHydraulics{FT}
    "[`AbstractReactionCenter`](@ref) type photosynthesis reaction center"
    PRC::Union{VJPReactionCenter{FT}, CytochromeReactionCenter{FT}}
    "[`AbstractPhotosynthesisModel`](@ref) type photosynthesis model"
    PSM::Union{C3VJPModel{FT}, C4VJPModel{FT}, C3CytochromeModel{FT}}
    "Stomatal model"
    SM::AbstractStomataModel{FT}
    "Leaf width"
    WIDTH::FT

    # prognostic variables that change with time
    "Stomatal conductance to water vapor `[mol m⁻² s⁻¹]`"
    g_H₂O_s::FT
    "Absorbed photosynthetically active radiation used for photosynthesis `[μmol m⁻² s⁻¹]`"
    ppar::FT
    "Current leaf temperature"
    t::FT

    # dignostic variables that change with time
    "Total leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_CO₂::FT
    "Boundary leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_CO₂_b::FT
    "Leaf internal CO₂ partial pressure `[Pa]`"
    p_CO₂_i::FT
    "Leaf surface CO₂ partial pressure `[Pa]`"
    p_CO₂_s::FT
    "Saturation H₂O vapor pressure, need to update with temperature and leaf water pressure `[Pa]`"
    p_H₂O_sat::FT

    # caches to speed up calculations
    "Last leaf temperature. If different from t, then make temperature correction"
    _t::FT
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
#     2022-Jun-15: add broadband as an option (default is false)
#     2022-Jun-29: add APAR_CAR as a field
#     2022-Jun-30: add SM as a field
#     2022-Jul-01: add G_LIMITS as a field
#
#######################################################################################################################################################################################################
"""

    Leaf{FT}(psm::String, wls::WaveLengthSet{FT} = WaveLengthSet{FT}(); broadband::Bool = false, colimit::Bool = false, ssm::Bool = true) where {FT<:AbstractFloat}

Constructor for `Leaf`, given
- `psm` Photosynthesis model type, must be `C3`, `C3Cytochrome`, or `C4`
- `wls` [`WaveLengthSet`](@ref) type structure that determines the dimensions of leaf parameters
- `broadband` Whether leaf biophysics is in broadband mode
- `colimit` Whether to colimit the photosynthetic rates and electron transport rates
- `ssm` Whether the flow rate is at steady state

---
# Examples
```julia
leaf_c3 = Leaf{Float64}("C3");
leaf_c4 = Leaf{Float64}("C4");
leaf_cy = Leaf{Float64}("C3Cytochrome");
leaf_c3 = Leaf{Float64}("C3"; colimit = true);
leaf_c4 = Leaf{Float64}("C4"; colimit = true);
leaf_cy = Leaf{Float64}("C3Cytochrome"; colimit = true);
wls = WaveLengthSet{FT}(collect(400:10:2500));
leaf_c3 = Leaf{Float64}("C3", wls);
leaf_c4 = Leaf{Float64}("C4", wls);
leaf_cy = Leaf{Float64}("C3Cytochrome", wls);
```
"""
Leaf{FT}(psm::String, wls::WaveLengthSet{FT} = WaveLengthSet{FT}(); broadband::Bool = false, colimit::Bool = false, ssm::Bool = true) where {FT<:AbstractFloat} = (
    @assert psm in ["C3", "C3Cytochrome", "C4"] "Photosynthesis model ID must be C3, C4, or C3Cytochrome!";

    if psm == "C3"
        _prc = VJPReactionCenter{FT}();
        _psm = C3VJPModel{FT}(colimit = colimit);
    elseif psm == "C3Cytochrome"
        _prc = CytochromeReactionCenter{FT}();
        _psm = C3CytochromeModel{FT}(colimit = colimit);
    elseif psm == "C4"
        _prc = VJPReactionCenter{FT}();
        _psm = C4VJPModel{FT}(colimit = colimit);
    end;

    if broadband
        _bio = BroadbandLeafBiophysics{FT}();
    else
        _bio = HyperspectralLeafBiophysics{FT}(wls);
    end;

    return Leaf{FT}(
                true,                               # APAR_CAR
                _bio,                               # BIO
                FT[0.01,0.3],                       # G_LIMITS
                LeafHydraulics{FT}(ssm = ssm),      # HS
                _prc,                               # PRC
                _psm,                               # PSM
                WangSM{FT}(),                       # SM
                FT(0.05),                           # WIDTH
                0.01,                               # g_H₂O_s
                1000,                               # ppar
                T_25(),                             # t
                0.01,                               # g_CO₂
                3.0,                                # g_CO₂_b
                20,                                 # p_CO₂_i
                40,                                 # p_CO₂_s
                saturation_vapor_pressure(T_25()),  # p_H₂O_sat
                0)                                  # _t
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
mutable struct Leaves1D{FT<:AbstractFloat}
    # parameters that do not change with time
    "[`BroadbandLeafBiophysics`](@ref) type leaf biophysical parameters"
    BIO::BroadbandLeafBiophysics{FT}
    "Minimal and maximum stomatal conductance for H₂O at 25 °C `[mol m⁻² s⁻¹]`"
    G_LIMITS::Vector{FT}
    "[`LeafHydraulics`](@ref) type leaf hydraulic system"
    HS::LeafHydraulics{FT}
    "[`LeafHydraulics`](@ref) type leaf hydraulic system used for other calculations (say sunlit and shaded leaf partitioning)"
    HS2::LeafHydraulics{FT}
    "[`AbstractReactionCenter`](@ref) type photosynthesis reaction center"
    PRC::Union{VJPReactionCenter{FT}, CytochromeReactionCenter{FT}}
    "[`AbstractPhotosynthesisModel`](@ref) type photosynthesis model"
    PSM::Union{C3VJPModel{FT}, C4VJPModel{FT}, C3CytochromeModel{FT}}
    "Stomatal model"
    SM::AbstractStomataModel{FT}
    "Leaf width"
    WIDTH::FT

    # prognostic variables that change with time
    "Stomatal conductance to water vapor `[mol m⁻² s⁻¹]`"
    g_H₂O_s::Vector{FT}
    "Absorbed photosynthetically active radiation used for photosynthesis `[μmol m⁻² s⁻¹]`"
    ppar::Vector{FT}
    "Current leaf temperature"
    t::Vector{FT}

    # dignostic variables that change with time
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_gross::Vector{FT}
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_net::Vector{FT}
    "Total leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_CO₂::Vector{FT}
    "Boundary leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_CO₂_b::Vector{FT}
    "Leaf internal CO₂ partial pressure `[Pa]`"
    p_CO₂_i::Vector{FT}
    "Leaf surface CO₂ partial pressure `[Pa]`"
    p_CO₂_s::Vector{FT}
    "Saturation H₂O vapor pressure, need to update with temperature and leaf water pressure `[Pa]`"
    p_H₂O_sat::FT
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-27: add constructor for Leaves1D
#     2022-Jun-27: make BIO BroadbandLeafBiophysics only
#     2022-Jun-28: add a_gross and a_net, make t a Vector, remove _t
#     2022-Jun-30: add a second HS2 for shaded leaves
#     2022-Jun-30: add SM as a field
#     2022-Jul-01: add G_LIMITS as a field
#
#######################################################################################################################################################################################################
"""

    Leaves1D{FT}(psm::String; colimit::Bool = false, ssm::Bool = true) where {FT<:AbstractFloat}

Constructor for `Leaves1D`, given
- `psm` Photosynthesis model type, must be `C3`, `C3Cytochrome`, or `C4`
- `colimit` Whether to colimit the photosynthetic rates and electron transport rates
- `ssm` Whether the flow rate is at steady state

---
# Examples
```julia
leaves_c3 = Leaves1D{Float64}("C3");
leaves_c4 = Leaves1D{Float64}("C4");
leaves_cy = Leaves1D{Float64}("C3Cytochrome");
leaves_c3 = Leaves1D{Float64}("C3"; colimit = true);
leaves_c4 = Leaves1D{Float64}("C4"; colimit = true);
leaves_cy = Leaves1D{Float64}("C3Cytochrome"; colimit = true);
```
"""
Leaves1D{FT}(psm::String; colimit::Bool = false, ssm::Bool = true) where {FT<:AbstractFloat} = (
    @assert psm in ["C3", "C3Cytochrome", "C4"] "Photosynthesis model ID must be C3, C4, or C3Cytochrome!";

    if psm == "C3"
        _prc = VJPReactionCenter{FT}();
        _psm = C3VJPModel{FT}(colimit = colimit);
    elseif psm == "C3Cytochrome"
        _prc = CytochromeReactionCenter{FT}();
        _psm = C3CytochromeModel{FT}(colimit = colimit);
    elseif psm == "C4"
        _prc = VJPReactionCenter{FT}();
        _psm = C4VJPModel{FT}(colimit = colimit);
    end;

    _bio = BroadbandLeafBiophysics{FT}();

    return Leaves1D{FT}(
                _bio,                               # BIO
                FT[0.01,0.3],                       # G_LIMITS
                LeafHydraulics{FT}(ssm = ssm),      # HS
                LeafHydraulics{FT}(ssm = ssm),      # HS2
                _prc,                               # PRC
                _psm,                               # PSM
                WangSM{FT}(),                       # SM
                FT(0.05),                           # WIDTH
                FT[0.01, 0.01],                     # g_H₂O_s
                FT[1000, 200],                      # ppar
                FT[T_25(), T_25()],                 # t
                zeros(FT,2),                        # a_gross
                zeros(FT,2),                        # a_net
                FT[0.01, 0.01],                     # g_CO₂
                FT[3.0, 3.0],                       # g_CO₂_b
                FT[20, 20],                         # p_CO₂_i
                FT[40, 40],                         # p_CO₂_s
                saturation_vapor_pressure(T_25()),  # p_H₂O_sat
    )
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
mutable struct Leaves2D{FT<:AbstractFloat}
    # parameters that do not change with time
    "Whether APAR absorbed by carotenoid is counted as PPAR"
    APAR_CAR::Bool
    "[`HyperspectralLeafBiophysics`](@ref) type leaf biophysical parameters"
    BIO::HyperspectralLeafBiophysics{FT}
    "Minimal and maximum stomatal conductance for H₂O at 25 °C `[mol m⁻² s⁻¹]`"
    G_LIMITS::Vector{FT}
    "[`LeafHydraulics`](@ref) type leaf hydraulic system"
    HS::LeafHydraulics{FT}
    "[`AbstractReactionCenter`](@ref) type photosynthesis reaction center"
    PRC::Union{VJPReactionCenter{FT}, CytochromeReactionCenter{FT}}
    "[`AbstractPhotosynthesisModel`](@ref) type photosynthesis model"
    PSM::Union{C3VJPModel{FT}, C4VJPModel{FT}, C3CytochromeModel{FT}}
    "Stomatal model"
    SM::AbstractStomataModel{FT}
    "Leaf width"
    WIDTH::FT

    # prognostic variables that change with time
    "Stomatal conductance to water vapor for shaded leaves `[mol m⁻² s⁻¹]`"
    g_H₂O_s_shaded::FT
    "Stomatal conductance to water vapor for sunlit leaves `[mol m⁻² s⁻¹]`"
    g_H₂O_s_sunlit::Matrix{FT}
    "Absorbed photosynthetically active radiation used for photosynthesis for shaded leaves `[μmol m⁻² s⁻¹]`"
    ppar_shaded::FT
    "Absorbed photosynthetically active radiation used for photosynthesis for sunlit leaves `[μmol m⁻² s⁻¹]`"
    ppar_sunlit::Matrix{FT}
    "Current leaf temperature"
    t::FT

    # dignostic variables that change with time
    "Gross photosynthetic rate for shaded leaves `[μmol m⁻² s⁻¹]`"
    a_gross_shaded::FT
    "Gross photosynthetic rate for sunlit leaves `[μmol m⁻² s⁻¹]`"
    a_gross_sunlit::Matrix{FT}
    "Net photosynthetic rate for shaded leaves `[μmol m⁻² s⁻¹]`"
    a_net_shaded::FT
    "Net photosynthetic rate for sunlit leaves `[μmol m⁻² s⁻¹]`"
    a_net_sunlit::Matrix{FT}
    "Total leaf diffusive conductance to CO₂ for shaded leaves `[mol m⁻² s⁻¹]`"
    g_CO₂_shaded::FT
    "Total leaf diffusive conductance to CO₂ for sunlit leaves `[mol m⁻² s⁻¹]`"
    g_CO₂_sunlit::Matrix{FT}
    "Boundary leaf diffusive conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_CO₂_b::FT
    "Leaf internal CO₂ partial pressure for shaded leaves `[Pa]`"
    p_CO₂_i_shaded::FT
    "Leaf internal CO₂ partial pressure for sunlit leaves `[Pa]`"
    p_CO₂_i_sunlit::Matrix{FT}
    "Leaf surface CO₂ partial pressure for shaded leaves `[Pa]`"
    p_CO₂_s_shaded::FT
    "Leaf surface CO₂ partial pressure for sunlit leaves `[Pa]`"
    p_CO₂_s_sunlit::Matrix{FT}
    "Saturation H₂O vapor pressure, need to update with temperature and leaf water pressure `[Pa]`"
    p_H₂O_sat::FT
    "Fluorescence quantum yield for shaded leaves `[-]`"
    ϕ_f_shaded::FT
    "Fluorescence quantum yield for sunlit leaves `[-]`"
    ϕ_f_sunlit::Matrix{FT}

    # caches to speed up calculations
    "Last leaf temperature. If different from t, then make temperature correction"
    _t::FT
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-27: add constructor for Leaves2D
#     2022-Jun-27: make BIO HyperspectralLeafBiophysics only
#     2022-Jun-27: add sunlit and shaded ppar to struct (remove the ppar in canopy radiation)
#     2022-Jun-28: add a_gross, a_net, and ϕ_f for sunlit and shaded leaves
#     2022-Jun-29: add APAR_CAR as a field
#     2022-Jun-30: add SM as a field
#     2022-Jul-01: add G_LIMITS as a field
#
#######################################################################################################################################################################################################
"""

    Leaves2D{FT}(psm::String, wls::WaveLengthSet{FT} = WaveLengthSet{FT}(); colimit::Bool = false, n_azi::Int = 36, n_incl::Int = 9, ssm::Bool = true) where {FT<:AbstractFloat}

Constructor for `Leaves2D`, given
- `psm` Photosynthesis model type, must be `C3`, `C3Cytochrome`, or `C4`
- `wls` [`WaveLengthSet`](@ref) type structure that determines the dimensions of leaf parameters
- `colimit` Whether to colimit the photosynthetic rates and electron transport rates
- `n_azi` Number of azimuth angles
- `n_incl` Number of inclination angles
- `ssm` Whether the flow rate is at steady state

---
# Examples
```julia
leaves_c3 = Leaves2D{Float64}("C3");
leaves_c4 = Leaves2D{Float64}("C4");
leaves_cy = Leaves2D{Float64}("C3Cytochrome");
leaves_c3 = Leaves2D{Float64}("C3"; colimit = true);
leaves_c4 = Leaves2D{Float64}("C4"; colimit = true);
leaves_cy = Leaves2D{Float64}("C3Cytochrome"; colimit = true);
wls = WaveLengthSet{FT}(collect(400:10:2500));
leaves_c3 = Leaves2D{Float64}("C3", wls);
leaves_c4 = Leaves2D{Float64}("C4", wls);
leaves_cy = Leaves2D{Float64}("C3Cytochrome", wls);
```
"""
Leaves2D{FT}(psm::String, wls::WaveLengthSet{FT} = WaveLengthSet{FT}(); colimit::Bool = false, n_azi::Int = 36, n_incl::Int = 9, ssm::Bool = true) where {FT<:AbstractFloat} = (
    @assert psm in ["C3", "C3Cytochrome", "C4"] "Photosynthesis model ID must be C3, C4, or C3Cytochrome!";

    if psm == "C3"
        _prc = VJPReactionCenter{FT}();
        _psm = C3VJPModel{FT}(colimit = colimit);
    elseif psm == "C3Cytochrome"
        _prc = CytochromeReactionCenter{FT}();
        _psm = C3CytochromeModel{FT}(colimit = colimit);
    elseif psm == "C4"
        _prc = VJPReactionCenter{FT}();
        _psm = C4VJPModel{FT}(colimit = colimit);
    end;

    _bio = HyperspectralLeafBiophysics{FT}(wls);

    return Leaves2D{FT}(
                true,                               # APAR_CAR
                _bio,                               # BIO
                FT[0.01,0.3],                       # G_LIMITS
                LeafHydraulics{FT}(ssm = ssm),      # HS
                _prc,                               # PRC
                _psm,                               # PSM
                WangSM{FT}(),                       # SM
                FT(0.05),                           # WIDTH
                0.01,                               # g_H₂O_s_shaded
                zeros(FT,n_incl,n_azi) .* FT(0.01), # g_H₂O_s_sunlit
                100,                                # ppar_shaded
                zeros(FT,n_incl,n_azi) .* 100,      # ppar_sunlit
                T_25(),                             # t
                0,                                  # a_gross_shaded
                zeros(FT,n_incl,n_azi),             # a_gross_sunlit
                0,                                  # a_net_shaded
                zeros(FT,n_incl,n_azi),             # a_net_sunlit
                0.01,                               # g_CO₂_shaded
                zeros(FT,n_incl,n_azi) .* FT(0.01), # g_CO₂_sunlit
                3.0,                                # g_CO₂_b
                20,                                 # p_CO₂_i_shaded
                zeros(FT,n_incl,n_azi) .+ 20,       # p_CO₂_i_sunlit
                40,                                 # p_CO₂_s_shaded
                zeros(FT,n_incl,n_azi) .+ 40,       # p_CO₂_s_sunlit
                saturation_vapor_pressure(T_25()),  # p_H₂O_sat
                0,                                  # ϕ_f_shaded
                zeros(FT,n_incl,n_azi),             # ϕ_f_sunlit
                0)                                  # _t
);
