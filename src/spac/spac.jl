#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-May-25: add abstract type for soil-plant-air continuum
#     2022-Jun-29: rename grass, palm, and tree SPAC to ML*SPAC
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractSPACSystem:
- [`MonoElementSPAC`](@ref)
- [`MonoMLGrassSPAC`](@ref)
- [`MonoMLPalmSPAC`](@ref)
- [`MonoMLTreeSPAC`](@ref)

"""
abstract type AbstractSPACSystem{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-25: toy SPAC system
#     2022-May-25: use Root and Stem structures with temperatures
#     2022-Jun-29: add AirLayer to SPAC
#     2022-Jul-14: add Meteorology to SPAC
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for simplest SPAC system

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct MonoElementSPAC{FT<:AbstractFloat} <: AbstractSPACSystem{FT}
    # Embedded structures
    "Air conditions"
    AIR::AirLayer{FT} = AirLayer{FT}()
    "Leaf system"
    LEAF::Leaf{FT} = Leaf{FT}()
    "Meteorology information"
    METEO::Meteorology{FT} = Meteorology{FT}()
    "Root system"
    ROOT::Root{FT} = Root{FT}()
    "Soil component"
    SOIL::Soil{FT} = Soil{FT}(ZS = FT[0, -1])
    "Stem system"
    STEM::Stem{FT} = Stem{FT}()

    # Cache variables
    "Relative hydraulic conductance"
    _krs::Vector{FT} = ones(FT, 4)
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-25: SPAC system for monospecies grass
#     2022-May-25: use Root and Stem structures with temperatures
#     2022-May-31: rename _qs to _fs
#     2022-Jun-29: rename struct to MonoMLPalmTreeSPAC, and use Leaves2D
#     2022-Jun-29: add CANOPY, Z, AIR, WLSET, LHA, ANGLES, SOIL, RAD_LW, RAD_SW, Φ_PHOTON to SPAC
#     2022-Jul-14: add Meteorology to SPAC
#     2022-Aug-30: remove LHA and WLSET
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for monospecies grass SPAC system

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct MonoMLGrassSPAC{FT<:AbstractFloat} <: AbstractSPACSystem{FT}
    # Dimensions
    "Dimension of air layers"
    DIM_AIR::Int = 25
    "Dimension of canopy layers"
    DIM_LAYER::Int = 20
    "Dimension of root layers"
    DIM_ROOT::Int = 5

    # General model information
    "Whether to convert energy to photons when computing fluorescence"
    Φ_PHOTON::Bool = true

    # Geographical information
    "Elevation"
    ELEVATION::FT = 32.0
    "Latitude"
    LATITUDE::FT = 33.173
    "Longitude"
    LONGITUDE::FT = 115.4494

    # Geometry information
    "Corresponding air layer per canopy layer"
    LEAVES_INDEX::Vector{Int} = collect(Int, 1:DIM_LAYER)
    "Corresponding soil layer per root layer"
    ROOTS_INDEX::Vector{Int} = collect(Int, 1:DIM_ROOT)
    "Depth and height information `[m]`"
    Z::Vector{FT} = FT[-0.2, 0.5]
    "Air boundaries `[m]`"
    Z_AIR::Vector{FT} = collect(FT, 0:0.05:1)

    # Embedded structures
    "Air for each layer (may be more than canopy layer)" # TODO: add air ZS
    AIR::Vector{AirLayer{FT}} = AirLayer{FT}[AirLayer{FT}() for _i in 1:DIM_AIR]
    "Sun sensor geometry"
    ANGLES::SunSensorGeometry{FT} = SunSensorGeometry{FT}()
    "Canopy used for radiation calculations"
    CANOPY::HyperspectralMLCanopy{FT} = HyperspectralMLCanopy{FT}()
    "Leaf per layer"
    LEAVES::Vector{Leaves2D{FT}} = Leaves2D{FT}[Leaves2D{FT}() for _i in 1:DIM_LAYER]
    "Meteorology information"
    METEO::Meteorology{FT} = Meteorology{FT}()
    "Downwelling longwave radiation `[W m⁻²]`"
    RAD_LW::FT = 100
    "Downwelling shortwave radiation"
    RAD_SW::HyperspectralRadiation{FT} = HyperspectralRadiation{FT}()
    "Root hydraulic system"
    ROOTS::Vector{Root{FT}} = Root{FT}[Root{FT}() for _i in 1:DIM_ROOT]
    "Soil component"
    SOIL::Soil{FT} = Soil{FT}()

    # Cache variables
    "Flow rate per root layer"
    _fs::Vector{FT} = zeros(FT, DIM_ROOT)
    "Conductances for each root layer at given flow"
    _ks::Vector{FT} = zeros(FT, DIM_ROOT)
    "Pressure for each root layer at given flow"
    _ps::Vector{FT} = zeros(FT, DIM_ROOT)
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-25: SPAC system for monospecies palm
#     2022-May-25: use Root and Stem structures with temperatures
#     2022-May-31: rename _qs to _fs
#     2022-Jun-29: rename struct to MonoMLPalmTreeSPAC, and use Leaves2D
#     2022-Jun-29: add CANOPY, Z, AIR, WLSET, LHA, ANGLES, SOIL, RAD_LW, RAD_SW, Φ_PHOTON to SPAC
#     2022-Jul-14: add Meteorology to SPAC
#     2022-Aug-30: remove LHA and WLSET
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for monospecies palm SPAC system (with trunk)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct MonoMLPalmSPAC{FT<:AbstractFloat} <: AbstractSPACSystem{FT}
    # Dimensions
    "Dimension of air layers"
    DIM_AIR::Int = 25
    "Dimension of canopy layers"
    DIM_LAYER::Int = 20
    "Dimension of root layers"
    DIM_ROOT::Int = 5

    # General model information
    "Whether to convert energy to photons when computing fluorescence"
    Φ_PHOTON::Bool = true

    # Geographical information
    "Elevation"
    ELEVATION::FT = 32.0
    "Latitude"
    LATITUDE::FT = 33.173
    "Longitude"
    LONGITUDE::FT = 115.4494

    # Geometry information
    "Corresponding air layer per canopy layer"
    LEAVES_INDEX::Vector{Int} = collect(Int, 1:DIM_LAYER)
    "Corresponding soil layer per root layer"
    ROOTS_INDEX::Vector{Int} = collect(Int, 1:DIM_ROOT)
    "Depth and height information `[m]`"
    Z::Vector{FT} = FT[-1, 6, 12]
    "Air boundaries `[m]`"
    Z_AIR::Vector{FT} = collect(FT, 0:0.2:13)

    # Embedded structures
    "Air for each layer (more than canopy layer)"
    AIR::Vector{AirLayer{FT}} = AirLayer{FT}[AirLayer{FT}() for _i in 1:DIM_AIR]
    "Sun sensor geometry"
    ANGLES::SunSensorGeometry{FT} = SunSensorGeometry{FT}()
    "Canopy used for radiation calculations"
    CANOPY::HyperspectralMLCanopy{FT} = HyperspectralMLCanopy{FT}()
    "Leaf per layer"
    LEAVES::Vector{Leaves2D{FT}} = Leaves2D{FT}[Leaves2D{FT}() for _i in 1:DIM_LAYER]
    "Meteorology information"
    METEO::Meteorology{FT} = Meteorology{FT}()
    "Downwelling longwave radiation `[W m⁻²]`"
    RAD_LW::FT = 100
    "Downwelling shortwave radiation"
    RAD_SW::HyperspectralRadiation{FT} = HyperspectralRadiation{FT}()
    "Root hydraulic system"
    ROOTS::Vector{Root{FT}} = Root{FT}[Root{FT}() for _i in 1:DIM_ROOT]
    "Soil component"
    SOIL::Soil{FT} = Soil{FT}()
    "Trunk hydraulic system"
    TRUNK::Stem{FT} = Stem{FT}()

    # Cache variables
    "Flow rate per root layer"
    _fs::Vector{FT} = zeros(FT, DIM_ROOT)
    "Conductances for each root layer at given flow"
    _ks::Vector{FT} = zeros(FT, DIM_ROOT)
    "Pressure for each root layer at given flow"
    _ps::Vector{FT} = zeros(FT, DIM_ROOT)
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-25: SPAC system for monospecies tree
#     2022-May-25: use Root and Stem structures with temperatures
#     2022-May-31: rename _qs to _fs
#     2022-Jun-29: rename struct to MonoMLTreeSPAC, and use Leaves2D
#     2022-Jun-29: add CANOPY, Z, AIR, WLSET, LHA, ANGLES, SOIL, RAD_LW, RAD_SW, Φ_PHOTON to SPAC
#     2022-Jul-14: add Meteorology to SPAC
#     2022-Aug-30: remove LHA and WLSET
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for monospecies tree SPAC system (with trunk and branches)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct MonoMLTreeSPAC{FT<:AbstractFloat} <: AbstractSPACSystem{FT}
    # dimensions
    "Dimension of air layers"
    DIM_AIR::Int = 25
    "Dimension of canopy layers"
    DIM_LAYER::Int = 20
    "Dimension of root layers"
    DIM_ROOT::Int = 5

    # General model information
    "Whether to convert energy to photons when computing fluorescence"
    Φ_PHOTON::Bool = true

    # Geometry information
    "Corresponding air layer per canopy layer"
    LEAVES_INDEX::Vector{Int} = collect(Int, 1:DIM_LAYER)
    "Corresponding soil layer per root layer"
    ROOTS_INDEX::Vector{Int} = collect(Int, 1:DIM_ROOT)
    "Depth and height information `[m]`"
    Z::Vector{FT} = FT[-1, 6, 12]
    "Air boundaries `[m]`"
    Z_AIR::Vector{FT} = collect(FT, 0:0.2:13)

    # Geographical information
    "Elevation"
    ELEVATION::FT = 32.0
    "Latitude"
    LATITUDE::FT = 33.173
    "Longitude"
    LONGITUDE::FT = 115.4494

    # Embedded structures
    "Air for each layer (more than canopy layer)"
    AIR::Vector{AirLayer{FT}} = AirLayer{FT}[AirLayer{FT}() for _i in 1:DIM_AIR]
    "Sun sensor geometry"
    ANGLES::SunSensorGeometry{FT} = SunSensorGeometry{FT}()
    "Branch hydraulic system"
    BRANCHES::Vector{Stem{FT}} = Stem{FT}[Stem{FT}() for _i in 1:DIM_LAYER]
    "Canopy used for radiation calculations"
    CANOPY::HyperspectralMLCanopy{FT} = HyperspectralMLCanopy{FT}()
    "Leaf per layer"
    LEAVES::Vector{Leaves2D{FT}} = Leaves2D{FT}[Leaves2D{FT}() for _i in 1:DIM_LAYER]
    "Meteorology information"
    METEO::Meteorology{FT} = Meteorology{FT}()
    "Downwelling longwave radiation `[W m⁻²]`"
    RAD_LW::FT = 100
    "Downwelling shortwave radiation"
    RAD_SW::HyperspectralRadiation{FT} = HyperspectralRadiation{FT}()
    "Root hydraulic system"
    ROOTS::Vector{Root{FT}} = Root{FT}[Root{FT}() for _i in 1:DIM_ROOT]
    "Soil component"
    SOIL::Soil{FT} = Soil{FT}()
    "Trunk hydraulic system"
    TRUNK::Stem{FT} = Stem{FT}()

    # Cache variables
    "Flow rate per root layer"
    _fs::Vector{FT} = zeros(FT, DIM_ROOT)
    "Conductances for each root layer at given flow"
    _ks::Vector{FT} = zeros(FT, DIM_ROOT)
    "Pressure for each root layer at given flow"
    _ps::Vector{FT} = zeros(FT, DIM_ROOT)
end


#=
TODO: move it to SoilPlantAirContinuum.jl
#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-May-25: add constructor function
#     2022-May-25: use Root and Stem structures with temperatures
#     2022-May-31: rename _qs to _fs
#     2022-May-31: add steady state mode option to input options
#     2022-Jun-15: fix documentation
#     2022-Jun-29: rename struct to MonoMLPalmTreeSPAC, and use Leaves2D
#     2022-Jun-29: add CANOPY, Z, AIR, WLSET, LHA, ANGLES, SOIL, RAD_LW, RAD_SW, Φ_PHOTON to SPAC
#     2022-Jul-14: add area to constructor function
#
#######################################################################################################################################################################################################
"""

    MonoMLGrassSPAC{FT}(
                psm::String,
                area::Number = 100,
                wls::WaveLengthSet{FT} = WaveLengthSet{FT}();
                zs::Vector = [-0.2,0.5],
                zss::Vector = collect(0:-0.1:-1),
                zas::Vector = collect(0:0.05:1),
                ssm::Bool = true
    ) where {FT<:AbstractFloat}

Construct a SPAC system for monospecies grass system, given
- `psm` Photosynthesis model, must be C3, C4, or C3Cytochrome
- `wls` [`WaveLengthSet`](@ref) type structure that determines the dimensions of leaf parameters
- `zs` Vector of Maximal root depth (negative value), and canopy height
- `zss` Vector of soil layer boundaries starting from 0
- `zas` Vector of air layer boundaries starting from 0
- `broadband` Whether leaf biophysics is in broadband mode
- `ssm` Whether the flow rate is at steady state

---
# Examples
```julia
spac = MonoMLGrassSPAC{Float64}("C3");
```
"""
MonoMLGrassSPAC{FT}(
            psm::String,
            area::Number = 100,
            wls::WaveLengthSet{FT} = WaveLengthSet{FT}();
            zs::Vector = [-0.2,0.5],
            zss::Vector = collect(0:-0.1:-1),
            zas::Vector = collect(0:0.05:1),
            ssm::Bool = true
) where {FT<:AbstractFloat} = (
    @assert psm in ["C3", "C4", "C3Cytochrome"] "Photosynthesis model must be within [C3, C4, C3CytochromeModel]";

    # determine how many layers of roots
    _n_root = 0;
    _r_inds = Int[];
    for _i in eachindex(zss)
        if zss[_i] > zs[1]
            _n_root += 1;
            push!(_r_inds, _i);
        else
            break
        end;
    end;

    # determine how many layers of canopy
    _n_canopy = 0;
    _c_inds = Int[];
    for _i in eachindex(zas)
        if zas[_i] < zs[2]
            _n_canopy += 1;
            push!(_c_inds, _i);
        else
            break
        end;
    end;

    # create evenly distributed root system for now
    _roots = Root{FT}[];
    for _i in _r_inds
        _Δh = abs(max(zss[_i+1], zs[1]) + zss[_i]) / 2;
        _hs = RootHydraulics{FT}(AREA = 1/_n_root, K_X = 25/_n_root, ΔH = _Δh);
        if !ssm _hs.FLOW = NonSteadyStateFlow{FT}(DIM_CAPACITY = 5) end;
        _rt = Root{FT}();
        _rt.HS = _hs;
        push!(_roots, _rt);
    end;

    # create leaves from bottom canopy (trunk) to top canopy
    _leaves = [Leaves2D{FT}(psm, wls) for _i in 1:_n_canopy];
    for _leaf in _leaves
        _leaf.HS.AREA = 1500 / _n_canopy;
    end;

    # create canopy to use for radiative transfer
    _canopy = HyperspectralMLCanopy{FT}();

    # create air layers for all provided layers from bottom to top
    _airs = [AirLayer{FT}() for _i in 1:length(zas)-1];

    # create leaf hyperspectral absorption features
    _lha = HyperspectralAbsorption{FT}();

    # create sun sensor geometry
    _angles = SunSensorGeometry{FT}();

    # create soil
    _soil = Soil{FT}(ZS = zss);

    # create shortwave radiation
    _rad_sw = HyperspectralRadiation{FT}();

    # return plant
    return MonoMLGrassSPAC{FT}(
                _airs,              # AIR
                _angles,            # ANGLES
                _canopy,            # CANOPY
                _leaves,            # LEAVES
                _c_inds,            # LEAVES_INDEX
                _lha,               # LHA
                Meteorology{FT}(),  # METEO
                _n_canopy,          # N_CANOPY
                _n_root,            # N_ROOT
                100,                # RAD_LW
                _rad_sw,            # RAD_SW
                _roots,             # ROOTS
                _r_inds,            # ROOTS_INDEX
                _soil,              # SOIL
                wls,                # WLSET
                FT.(zs),            # Z
                FT.(zas),           # Z_AIR
                true,               # Φ_PHOTON
                zeros(FT,_n_root),  # _fs
                zeros(FT,_n_root),  # _ks
                zeros(FT,_n_root)   # _ps
    )
);


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-May-25: add constructor function
#     2022-May-25: use Root and Stem structures with temperatures
#     2022-May-31: rename _qs to _fs
#     2022-May-31: add steady state mode option to input options
#     2022-Jun-15: fix documentation
#     2022-Jun-29: rename struct to MonoMLPalmTreeSPAC, and use Leaves2D
#     2022-Jun-29: add CANOPY, Z, AIR, WLSET, LHA, ANGLES, SOIL, RAD_LW, RAD_SW, Φ_PHOTON to SPAC
#     2022-Jul-14: add area to constructor function
#
#######################################################################################################################################################################################################
"""

    MonoMLPalmSPAC{FT}(
                psm::String,
                area::Number = 100,
                wls::WaveLengthSet{FT} = WaveLengthSet{FT}();
                zs::Vector = [-1,6,12],
                zss::Vector = collect(0:-0.25:-2),
                zas::Vector = collect(0:0.2:13),
                ssm::Bool = true
    ) where {FT<:AbstractFloat}

Construct a SPAC system for monospecies palm system, given
- `psm` Photosynthesis model, must be C3 or C3Cytochrome
- `wls` [`WaveLengthSet`](@ref) type structure that determines the dimensions of leaf parameters
- `zs` Vector of Maximal root depth (negative value), trunk height, and canopy height
- `zss` Vector of soil layer boundaries starting from 0
- `zas` Vector of air layer boundaries starting from 0
- `ssm` Whether the flow rate is at steady state

---
# Examples
```julia
spac = MonoMLPalmSPAC{Float64}("C3");
```
"""
MonoMLPalmSPAC{FT}(
            psm::String,
            area::Number = 100,
            wls::WaveLengthSet{FT} = WaveLengthSet{FT}();
            zs::Vector = [-1,6,12],
            zss::Vector = collect(0:-0.25:-2),
            zas::Vector = collect(0:0.2:13),
            ssm::Bool = true
) where {FT<:AbstractFloat} = (
    @assert psm in ["C3", "C4", "C3Cytochrome"] "Photosynthesis model must be within [C3, C4, C3CytochromeModel]";

    # determine how many layers of roots
    _n_root = 0;
    _r_inds = Int[];
    for _i in eachindex(zss)
        if zss[_i] > zs[1]
            _n_root += 1;
            push!(_r_inds, _i);
        else
            break
        end;
    end;

    # determine how many layers of canopy
    _n_canopy= 0;
    _c_inds = Int[];
    for _i in 1:(length(zas)-1)
        # if the entire canopy is within the same layer
        if zas[_i] <= zs[2] < zs[3] <= zas[_i+1]
            _n_canopy += 1;
            push!(_c_inds, _i);
            break

        # if the lower canopy boundary (trunk top) is within the layer
        elseif zas[_i] < zs[2] < zas[_i+1] < zs[3]
            _n_canopy += 1;
            push!(_c_inds, _i);

        # if the layer is within the canopy
        elseif zs[2] <= zas[_i] < zas[_i+1] < zs[3]
            _n_canopy += 1;
            push!(_c_inds, _i);

        # if the upper canopy boundary is within the layer
        elseif zs[2] <= zas[_i] <= zs[3] < zas[_i+1]
            _n_canopy += 1;
            push!(_c_inds, _i-1);
            break

        # if the entire canopy is below the layer
        elseif zs[3] <= zas[_i] < zas[_i+1]
            break
        end;
    end;

    # create evenly distributed root system from top soil to deep soil
    _roots = Root{FT}[];
    for _i in _r_inds
        _Δh = abs(max(zss[_i+1], zs[1]) + zss[_i]) / 2;
        _hs = RootHydraulics{FT}(AREA = 1/_n_root, K_X = 25/_n_root, ΔH = _Δh);
        if !ssm _hs.FLOW = NonSteadyStateFlow{FT}(DIM_CAPACITY = 5) end;
        _rt = Root{FT}();
        _rt.HS = _hs;
        push!(_roots, _rt);
    end;

    # create trunk with the height of entire canopy
    _hs = StemHydraulics{FT}(L = zs[3], ΔH = zs[3]);
    if !ssm
        _hs.FLOW = NonSteadyStateFlow{FT}(DIM_CAPACITY = 5);
    end;
    _trunk = Stem{FT}();
    _trunk.HS = _hs;

    # create leaves from bottom canopy (trunk) to top canopy
    _leaves = [Leaves2D{FT}(psm, wls) for _i in 1:_n_canopy];
    for _leaf in _leaves
        _leaf.HS.AREA = 1500 / _n_canopy;
    end;

    # create canopy to use for radiative transfer
    _canopy = HyperspectralMLCanopy{FT}();

    # create air layers for all provided layers from bottom to top
    _airs = [AirLayer{FT}() for _i in 1:length(zas)-1];

    # create leaf hyperspectral absorption features
    _lha = HyperspectralAbsorption{FT}();

    # create sun sensor geometry
    _angles = SunSensorGeometry{FT}();

    # create soil
    _soil = Soil{FT}(ZS = zss);

    # create shortwave radiation
    _rad_sw = HyperspectralRadiation{FT}();

    # return plant
    return MonoMLPalmSPAC{FT}(
                _airs,              # AIR
                _angles,            # ANGLES
                _canopy,            # CANOPY
                _leaves,            # LEAVES
                _c_inds,            # LEAVES_INDEX
                _lha,               # LHA
                Meteorology{FT}(),  # METEO
                _n_canopy,          # N_CANOPY
                _n_root,            # N_ROOT
                100,                # RAD_LW
                _rad_sw,            # RAD_SW
                _roots,             # ROOTS
                _r_inds,            # ROOTS_INDEX
                _soil,              # SOIL
                _trunk,             # TRUNK
                wls,                # WLSET
                FT.(zs),            # Z
                FT.(zas),           # Z_AIR
                true,               # Φ_PHOTON
                zeros(FT,_n_root),  # _fs
                zeros(FT,_n_root),  # _ks
                zeros(FT,_n_root)   # _ps
    )
);


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-May-25: add constructor function
#     2022-May-25: use Root and Stem structures with temperatures
#     2022-May-31: rename _qs to _fs
#     2022-May-31: add steady state mode option to input options
#     2022-Jun-15: fix documentation
#     2022-Jun-29: rename struct to MonoMLTreeSPAC, and use Leaves2D
#     2022-Jun-29: add CANOPY, Z, AIR, WLSET, LHA, ANGLES, SOIL, RAD_LW, RAD_SW, Φ_PHOTON to SPAC
#     2022-Jul-14: add area to constructor function
#
#######################################################################################################################################################################################################
"""

    MonoMLTreeSPAC{FT}(
                psm::String,
                area::Number = 100,
                wls::WaveLengthSet{FT} = WaveLengthSet{FT}();
                zs::Vector = [-1,6,12],
                zss::Vector = collect(0:-0.25:-2),
                zas::Vector = collect(0:0.5:13),
                ssm::Bool = true
    ) where {FT<:AbstractFloat}

Construct a SPAC system for monospecies tree system, given
- `psm` Photosynthesis model, must be C3, C4, or C3Cytochrome (note: there are C4 shrubs)
- `wls` [`WaveLengthSet`](@ref) type structure that determines the dimensions of leaf parameters
- `zs` Vector of Maximal root depth (negative value), trunk height, and canopy height
- `zss` Vector of soil layer boundaries starting from 0
- `zas` Vector of air layer boundaries starting from 0
- `ssm` Whether the flow rate is at steady state

---
# Examples
```julia
spac = MonoMLTreeSPAC{Float64}("C3");
```
"""
MonoMLTreeSPAC{FT}(
            psm::String,
            area::Number = 100,
            wls::WaveLengthSet{FT} = WaveLengthSet{FT}();
            zs::Vector = [-1,6,12],
            zss::Vector = collect(0:-0.25:-2),
            zas::Vector = collect(0:0.5:13),
            ssm::Bool = true
) where {FT<:AbstractFloat} = (
    @assert psm in ["C3", "C4", "C3Cytochrome"] "Photosynthesis model must be within [C3, C4, C3CytochromeModel]";

    # determine how many layers of roots
    _n_root = 0;
    _r_inds = Int[];
    for _i in eachindex(zss)
        if zss[_i] > zs[1]
            _n_root += 1;
            push!(_r_inds, _i);
        else
            break
        end;
    end;

    # determine how many layers of canopy
    _n_canopy= 0;
    _c_inds = Int[];
    for _i in 1:(length(zas)-1)
        # if the entire canopy is within the same layer
        if zas[_i] <= zs[2] < zs[3] <= zas[_i+1]
            _n_canopy += 1;
            push!(_c_inds, _i);
            break

        # if the lower canopy boundary (trunk top) is within the layer
        elseif zas[_i] < zs[2] < zas[_i+1] < zs[3]
            _n_canopy += 1;
            push!(_c_inds, _i);

        # if the layer is within the canopy
        elseif zs[2] <= zas[_i] < zas[_i+1] < zs[3]
            _n_canopy += 1;
            push!(_c_inds, _i);

        # if the upper canopy boundary is within the layer
        elseif zs[2] <= zas[_i] <= zs[3] < zas[_i+1]
            _n_canopy += 1;
            push!(_c_inds, _i-1);
            break

        # if the entire canopy is below the layer
        elseif zs[3] <= zas[_i] < zas[_i+1]
            break
        end;
    end;

    # create evenly distributed root system from top soil to deep soil
    _roots = Root{FT}[];
    for _i in _r_inds
        _Δh = abs(max(zss[_i+1], zs[1]) + zss[_i]) / 2;
        _hs = RootHydraulics{FT}(AREA = 1/_n_root, K_X = 25/_n_root, ΔH = _Δh);
        if !ssm _hs.FLOW = NonSteadyStateFlow{FT}(DIM_CAPACITY = 5) end;
        _rt = Root{FT}();
        _rt.HS = _hs;
        push!(_roots, _rt);
    end;

    # create trunk
    _hs = StemHydraulics{FT}(L = zs[2], ΔH = zs[2]);
    if !ssm
        _hs.FLOW = NonSteadyStateFlow{FT}(DIM_CAPACITY = 5);
    end;
    _trunk = Stem{FT}();
    _trunk.HS = _hs;

    # create branches from bottom canopy (trunk) to top canopy
    _branches = Stem{FT}[];
    for _i in _c_inds
        _Δh = (max(zas[_i], zs[2]) + min(zas[_i+1], zs[3])) / 2 - zs[2];
        _hs = StemHydraulics{FT}(AREA = 1/_n_canopy, ΔH = _Δh);
        if !ssm
            _hs.FLOW = NonSteadyStateFlow{FT}(DIM_CAPACITY = 5);
        end;
        _st = Stem{FT}();
        _st.HS = _hs;
        push!(_branches, _st);
    end;

    # create leaves from bottom canopy (trunk) to top canopy
    _leaves = [Leaves2D{FT}(psm, wls) for _i in 1:_n_canopy];
    for _leaf in _leaves
        _leaf.HS.AREA = 1500 / _n_canopy;
    end;

    # create canopy to use for radiative transfer
    _canopy = HyperspectralMLCanopy{FT}();

    # create air layers for all provided layers from bottom to top
    _airs = [AirLayer{FT}() for _i in 1:length(zas)-1];

    # create leaf hyperspectral absorption features
    _lha = HyperspectralAbsorption{FT}();

    # create sun sensor geometry
    _angles = SunSensorGeometry{FT}();

    # create soil
    _soil = Soil{FT}(ZS = zss);

    # create shortwave radiation
    _rad_sw = HyperspectralRadiation{FT}();

    # return plant
    return MonoMLTreeSPAC{FT}(
                _airs,              # AIR
                _angles,            # ANGLES
                _branches,          # BRANCHES
                _canopy,            # CANOPY
                _leaves,            # LEAVES
                _c_inds,            # LEAVES_INDEX
                _lha,               # LHA
                Meteorology{FT}(),  # METEO
                _n_canopy,          # N_CANOPY
                _n_root,            # N_ROOT
                100,                # RAD_LW
                _rad_sw,            # RAD_SW
                _roots,             # ROOTS
                _r_inds,            # ROOTS_INDEX
                _soil,              # SOIL
                _trunk,             # TRUNK
                wls,                # WLSET
                FT.(zs),            # Z
                FT.(zas),           # Z_AIR
                true,               # Φ_PHOTON
                zeros(FT,_n_root),  # _fs
                zeros(FT,_n_root),  # _ks
                zeros(FT,_n_root)   # _ps
    )
);
=#
