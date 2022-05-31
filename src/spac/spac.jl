#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-May-25: add abstract type for soil-plant-air continuum
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractSPACSystem:
- [`MonoElementSPAC`](@ref)
- [`MonoGrassSPAC`](@ref)
- [`MonoPalmSPAC`](@ref)
- [`MonoTreeSPAC`](@ref)
"""
abstract type AbstractSPACSystem{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-25: toy SPAC system
#     2022-May-25: use Root and Stem structures with temperatures
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for simplest SPAC system

# Fields

$(TYPEDFIELDS)

"""
mutable struct MonoElementSPAC{FT} <: AbstractSPACSystem{FT}
    # parameters that do not change with time
    "Leaf system"
    LEAF::Leaf{FT}
    "Root system"
    ROOT::Root{FT}
    "Stem system"
    STEM::Stem{FT}

    # caches to speed up calculations
    "Relative hydraulic conductance"
    _krs::Vector{FT}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-May-25: add constructor function
#     2022-May-25: add psm to constructor option
#     2022-May-25: use Root and Stem structures with temperatures
#     2022-May-31: add steady state mode option to input options
#
#######################################################################################################################################################################################################
"""

    MonoElementSPAC{FT}(psm::String; ssm::Bool = true) where {FT<:AbstractFloat}

Construct a `MonoElementSPAC` type toy SPAC system, given
- `psm` Photosynthesis model, must be C3, C4, or C3Cytochrome
- `ssm` Whether the flow rate is at steady state
"""
MonoElementSPAC{FT}(psm::String; ssm::Bool = true) where {FT<:AbstractFloat} = (
    @assert psm in ["C3", "C4", "C3Cytochrome"] "Photosynthesis model must be within [C3, C4, C3CytochromeModel]";

    return MonoElementSPAC{FT}(
                Leaf{FT}(psm; ssm = ssm),   # LEAF
                Root{FT}(ssm = ssm),        # ROOT
                Stem{FT}(ssm = ssm),        # STEM
                ones(FT,4)                  # _krs
    )
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-25: SPAC system for monospecies grass
#     2022-May-25: use Root and Stem structures with temperatures
#     2022-May-31: rename _qs to _fs
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for monospecies grass SPAC system

# Fields

$(TYPEDFIELDS)

"""
mutable struct MonoGrassSPAC{FT} <: AbstractSPACSystem{FT}
    # parameters that do not change with time
    "Leaf hydrualic system"
    LEAVES::Vector{Leaf{FT}}
    "Corresponding air layer per canopy layer"
    LEAVES_INDEX::Vector{Int}
    "Number of canopy layers"
    N_CANOPY::Int
    "Number of root layers"
    N_ROOT::Int
    "Root hydraulic system"
    ROOTS::Vector{Root{FT}}
    "Corresponding soil layer per root layer"
    ROOTS_INDEX::Vector{Int}

    # caches to speed up calculations
    "Flow rate per root layer"
    _fs::Vector{FT}
    "Conductances for each root layer at given flow"
    _ks::Vector{FT}
    "Pressure for each root layer at given flow"
    _ps::Vector{FT}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-May-25: add constructor function
#     2022-May-25: use Root and Stem structures with temperatures
#     2022-May-31: rename _qs to _fs
#     2022-May-31: add steady state mode option to input options
#
#######################################################################################################################################################################################################
"""

    MonoGrassSPAC{FT}(psm::String; zr::Number = -0.2, zc::Number = 0.5, zss::Vector = collect(0:-0.1:-1), zas::Vector = collect(0:0.05:1), ssm::Bool = true) where {FT<:AbstractFloat}

Construct a SPAC system for monospecies grass system, given
- `psm` Photosynthesis model, must be C3, C4, or C3Cytochrome
- `zr` Maximal root depth (negative value)
- `zc` Maximal canopy height (positive value)
- `zss` Vector of soil layer boundaries starting from 0
- `zas` Vector of air layer boundaries starting from 0
- `ssm` Whether the flow rate is at steady state

---
# Examples
```julia
spac = MonoGrassSPAC{Float64}();
spac = MonoGrassSPAC{Float64}(zr = -0.3, zc = 1, zss = collect(0:-0.1:-1), zas = collect(0:0.05:1.01));
```
"""
MonoGrassSPAC{FT}(psm::String; zr::Number = -0.2, zc::Number = 0.5, zss::Vector = collect(0:-0.1:-1), zas::Vector = collect(0:0.05:1), ssm::Bool = true) where {FT<:AbstractFloat} = (
    @assert psm in ["C3", "C4", "C3Cytochrome"] "Photosynthesis model must be within [C3, C4, C3CytochromeModel]";

    # determine how many layers of roots
    _n_root  = 0;
    _r_index = Int[];
    for i in eachindex(zss)
        _z = zss[i];
        if _z > zr
            _n_root += 1;
            push!(_r_index, i);
        else
            break
        end;
    end;

    # determine how many layers of canopy
    _n_canopy = 0;
    _c_index  = Int[];
    for i in eachindex(zas)
        _z = zas[i]
        if _z < zc
            _n_canopy += 1;
            push!(_c_index, i);
        elseif _z >= zc
            break
        end;
    end;

    # create evenly distributed root system for now
    _roots = Root{FT}[];
    for i in _r_index
        _Δh = abs(zss[i+1] + zss[i]) / 2;
        _rt = Root{FT}(RootHydraulics{FT}(area = 1/_n_root, k_x = 25/_n_root, Δh = _Δh, ssm = ssm), T_25());
        push!(_roots, _rt);
    end;

    # create leaves
    _leaves = [Leaf{FT}(psm; ssm = ssm) for i in 1:_n_canopy];
    for _leaf in _leaves
        _leaf.HS.AREA = 1500 / _n_canopy;
    end;

    # return plant
    return MonoGrassSPAC{FT}(
                _leaves,            # LEAVES
                _c_index,           # LEVES_INDEX
                _n_canopy,          # N_CANOPY
                _n_root,            # N_ROOT
                _roots,             # ROOTS
                _r_index,           # ROOTS_INDEX
                zeros(FT,_n_root),  # _fs
                zeros(FT,_n_root),  # _ks
                zeros(FT,_n_root)   # _ps
    )
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-25: SPAC system for monospecies palm
#     2022-May-25: use Root and Stem structures with temperatures
#     2022-May-31: rename _qs to _fs
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for monospecies palm SPAC system (with trunk)

# Fields

$(TYPEDFIELDS)

"""
mutable struct MonoPalmSPAC{FT} <: AbstractSPACSystem{FT}
    # parameters that do not change with time
    "Leaf hydrualic system"
    LEAVES::Vector{Leaf{FT}}
    "Corresponding air layer per canopy layer"
    LEAVES_INDEX::Vector{Int}
    "Number of canopy layers"
    N_CANOPY::Int
    "Number of root layers"
    N_ROOT::Int
    "Root hydraulic system"
    ROOTS::Vector{Root{FT}}
    "Corresponding soil layer per root layer"
    ROOTS_INDEX::Vector{Int}
    "Trunk hydraulic system"
    TRUNK::Stem{FT}

    # caches to speed up calculations
    "Flow rate per root layer"
    _fs::Vector{FT}
    "Conductances for each root layer at given flow"
    _ks::Vector{FT}
    "Pressure for each root layer at given flow"
    _ps::Vector{FT}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-May-25: add constructor function
#     2022-May-25: use Root and Stem structures with temperatures
#     2022-May-31: rename _qs to _fs
#     2022-May-31: add steady state mode option to input options
#
#######################################################################################################################################################################################################
"""

    MonoPalmSPAC{FT}(psm::String; zr::Number = -1, zt::Number = 10, zc::Number = 12, zss::Vector = collect(0:-0.25:-2), zas::Vector = collect(0:0.2:13), ssm::Bool = true) where {FT<:AbstractFloat}

Construct a SPAC system for monospecies palm system, given
- `psm` Photosynthesis model, must be C3 or C3Cytochrome
- `zr` Maximal root depth (negative value)
- `zc` Maximal canopy height (positive value)
- `zss` Vector of soil layer boundaries starting from 0
- `zas` Vector of air layer boundaries starting from 0
- `ssm` Whether the flow rate is at steady state

---
# Examples
```julia
spac = MonoPalmSPAC{Float64}();
spac = MonoPalmSPAC{Float64}(zr = -1, zt = 11, zc = 1, zss = collect(0:-0.1:-2), zas = collect(0:0.2:13));
```
"""
MonoPalmSPAC{FT}(psm::String; zr::Number = -1, zt::Number = 10, zc::Number = 12, zss::Vector = collect(0:-0.25:-2), zas::Vector = collect(0:0.2:13), ssm::Bool = true) where {FT<:AbstractFloat} = (
    @assert psm in ["C3", "C3Cytochrome"] "Photosynthesis model must be within [C3, C3CytochromeModel]";

    # determine how many layers of roots
    _n_root  = 0;
    _r_index = Int[];
    for i in eachindex(zss)
        _z = zss[i];
        if _z > zr
            _n_root += 1;
            push!(_r_index, i);
        else
            break
        end;
    end;

    # determine how many layers of canopy
    _n_canopy = 0;
    _c_index  = Int[];
    for i in 1:(length(zas)-1)
        _z0 = zas[i];
        _z1 = zas[i+1];
        if _z0 <= zt < zc <= _z1
            _n_canopy += 1;
            push!(_c_index, i);
            break
        elseif _z0 < zt < _z1 < zc
            _n_canopy += 1;
            push!(_c_index, i);
        elseif zt <= _z0 < _z1 < zc
            _n_canopy += 1;
            push!(_c_index, i);
        elseif zt <= _z0 <= zc < _z1
            _n_canopy += 1;
            push!(_c_index, i-1);
            break
        elseif zc <= _z0 < _z1
            break
        end;
    end;

    # create evenly distributed root system for now
    _roots = Root{FT}[];
    for i in _r_index
        _Δh = abs(zss[i+1] + zss[i]) / 2;
        _rt = Root{FT}(RootHydraulics{FT}(area = 1/_n_root, k_x = 25/_n_root, Δh = _Δh, ssm = ssm), T_25());
        push!(_roots, _rt);
    end;

    # create trunk
    _trunk = Stem{FT}(StemHydraulics{FT}(Δh = zt, Δl = zt, ssm = ssm), T_25());

    # create leaves
    _leaves = [Leaf{FT}(psm; ssm = ssm) for i in 1:_n_canopy];
    for _leaf in _leaves
        _leaf.HS.AREA = 1500 / _n_canopy;
    end;

    # return plant
    return MonoPalmSPAC{FT}(
                _leaves,            # LEAVES
                _c_index,           # LEVES_INDEX
                _n_canopy,          # N_CANOPY
                _n_root,            # N_ROOT
                _roots,             # ROOTS
                _r_index,           # ROOTS_INDEX
                _trunk,             # TRUNK
                zeros(FT,_n_root),  # _fs
                zeros(FT,_n_root),  # _ks
                zeros(FT,_n_root)   # _ps
    )
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-25: SPAC system for monospecies tree
#     2022-May-25: use Root and Stem structures with temperatures
#     2022-May-31: rename _qs to _fs
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for monospecies tree SPAC system (with trunk and branches)

# Fields

$(TYPEDFIELDS)

"""
mutable struct MonoTreeSPAC{FT} <: AbstractSPACSystem{FT}
    # parameters that do not change with time
    "Branch hydraulic system"
    BRANCHES::Vector{Stem{FT}}
    "Leaf hydrualic system"
    LEAVES::Vector{Leaf{FT}}
    "Corresponding air layer per canopy layer"
    LEAVES_INDEX::Vector{Int}
    "Number of canopy layers"
    N_CANOPY::Int
    "Number of root layers"
    N_ROOT::Int
    "Root hydraulic system"
    ROOTS::Vector{Root{FT}}
    "Corresponding soil layer per root layer"
    ROOTS_INDEX::Vector{Int}
    "Trunk hydraulic system"
    TRUNK::Stem{FT}

    # caches to speed up calculations
    "Flow rate per root layer"
    _fs::Vector{FT}
    "Conductances for each root layer at given flow"
    _ks::Vector{FT}
    "Pressure for each root layer at given flow"
    _ps::Vector{FT}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-May-25: add constructor function
#     2022-May-25: use Root and Stem structures with temperatures
#     2022-May-31: rename _qs to _fs
#     2022-May-31: add steady state mode option to input options
#
#######################################################################################################################################################################################################
"""

    MonoTreeSPAC{FT}(psm::String; zr::Number = -1, zt::Number = 10, zc::Number = 12, zss::Vector = collect(0:-0.25:-2), zas::Vector = collect(0:0.2:13), ssm::Bool = true) where {FT<:AbstractFloat}

Construct a SPAC system for monospecies tree system, given
- `psm` Photosynthesis model, must be C3, C4, or C3Cytochrome (note: there are C4 shrubs)
- `zr` Maximal root depth (negative value)
- `zc` Maximal canopy height (positive value)
- `zss` Vector of soil layer boundaries starting from 0
- `zas` Vector of air layer boundaries starting from 0
- `ssm` Whether the flow rate is at steady state

---
# Examples
```julia
spac = MonoTreeSPAC{Float64}();
spac = MonoTreeSPAC{Float64}(zr = -1, zt = 11, zc = 1, zss = collect(0:-0.1:-2), zas = collect(0:0.2:13));
```
"""
MonoTreeSPAC{FT}(psm::String; zr::Number = -1, zt::Number = 10, zc::Number = 12, zss::Vector = collect(0:-0.25:-2), zas::Vector = collect(0:0.2:13), ssm::Bool = true) where {FT<:AbstractFloat} = (
    @assert psm in ["C3", "C4", "C3Cytochrome"] "Photosynthesis model must be within [C3, C4, C3CytochromeModel]";

    # determine how many layers of roots
    _n_root  = 0;
    _r_index = Int[];
    for i in eachindex(zss)
        _z = zss[i];
        if _z > zr
            _n_root += 1;
            push!(_r_index, i);
        else
            break
        end;
    end;

    # determine how many layers of canopy
    _n_canopy = 0;
    _c_index  = Int[];
    for i in 1:(length(zas)-1)
        _z0 = zas[i];
        _z1 = zas[i+1];
        if _z0 <= zt < zc <= _z1
            _n_canopy += 1;
            push!(_c_index, i);
            break
        elseif _z0 < zt < _z1 < zc
            _n_canopy += 1;
            push!(_c_index, i);
        elseif zt <= _z0 < _z1 < zc
            _n_canopy += 1;
            push!(_c_index, i);
        elseif zt <= _z0 <= zc < _z1
            _n_canopy += 1;
            push!(_c_index, i-1);
            break
        elseif zc <= _z0 < _z1
            break
        end;
    end;

    # create evenly distributed root system for now
    _roots = Root{FT}[];
    for i in _r_index
        _Δh = abs(zss[i+1] + zss[i]) / 2;
        _rt = Root{FT}(RootHydraulics{FT}(area = 1/_n_root, k_x = 25/_n_root, Δh = _Δh, ssm = ssm), T_25());
        push!(_roots, _rt);
    end;

    # create trunk
    _trunk = Stem{FT}(StemHydraulics{FT}(Δh = zt, Δl = zt, ssm = ssm), T_25());

    # create branches
    _branches = Stem{FT}[];
    for i in _c_index
        _Δh = (zas[i] + max(zas[i+1], zt)) / 2 - zt;
        _st = Stem{FT}(StemHydraulics{FT}(area = 1/_n_canopy, Δh = _Δh, ssm = ssm), T_25());
        push!(_branches, _st);
    end;

    # create leaves
    _leaves = [Leaf{FT}(psm; ssm = ssm) for i in 1:_n_canopy];
    for _leaf in _leaves
        _leaf.HS.AREA = 1500 / _n_canopy;
    end;

    # return plant
    return MonoTreeSPAC{FT}(
                _branches,          # BRANCHES
                _leaves,            # LEAVES
                _c_index,           # LEVES_INDEX
                _n_canopy,          # N_CANOPY
                _n_root,            # N_ROOT
                _roots,             # ROOTS
                _r_index,           # ROOTS_INDEX
                _trunk,             # TRUNK
                zeros(FT,_n_root),  # _fs
                zeros(FT,_n_root),  # _ks
                zeros(FT,_n_root)   # _ps
    )
);
