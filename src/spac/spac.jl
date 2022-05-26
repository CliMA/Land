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
    "Leaf hydrualic system"
    LEAF::Leaf{FT}
    "Root hydraulic system"
    ROOT::RootHydraulics{FT}
    "Stem hydraulic system"
    STEM::StemHydraulics{FT}

    # caches to speed up calculations
    "Relative hydraulic conductance"
    _krs::Array{FT,1}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-May-25: add constructor function
#
#######################################################################################################################################################################################################
"""

    MonoElementSPAC{FT}() where {FT<:AbstractFloat}

Construct a `MonoElementSPAC` type toy SPAC system
"""
MonoElementSPAC{FT}() where {FT<:AbstractFloat} = (
    return MonoElementSPAC{FT}(
                Leaf{FT}(),             # LEAF
                RootHydraulics{FT}(),   # ROOT
                StemHydraulics{FT}(),   # STEM
                ones(FT,3)              # _krs
    )
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-25: SPAC system for monospecies grass
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
    ROOTS::Vector{RootHydraulics{FT}}
    "Corresponding soil layer per root layer"
    ROOTS_INDEX::Vector{Int}

    # caches to speed up calculations
    "Conductances for each root layer at given flow"
    _ks::Array{FT,1}
    "Pressure for each root layer at given flow"
    _ps::Array{FT,1}
    "Flow rate per root layer"
    _qs::Array{FT,1}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-May-25: add constructor function
#
#######################################################################################################################################################################################################
"""

    MonoGrassSPAC{FT}(psm::String; zr::Number = -0.2, zc::Number = 0.5, z_soil::Vector = collect(0:-0.1:-1), z_air::Vector = collect(0:0.05:1)) where {FT<:AbstractFloat}

Construct a SPAC system for monospecies grass system, given
- `psm` Photosynthesis model, must be C3, C4, or C3Cytochrome
- `zr` Maximal root depth (negative value)
- `zc` Maximal canopy height (positive value)
- `z_soil` Array of soil layer boundaries starting from 0
- `z_air` Array of air layer boundaries starting from 0

---
# Examples
```julia
spac = MonoGrassSPAC{Float64}();
spac = MonoGrassSPAC{Float64}(zr = -0.3, zc = 1, z_soil = collect(0:-0.1:-1), z_air = collect(0:0.05:1.01));
```
"""
MonoGrassSPAC{FT}(psm::String; zr::Number = -0.2, zc::Number = 0.5, z_soil::Vector = collect(0:-0.1:-1), z_air::Vector = collect(0:0.05:1)) where {FT<:AbstractFloat} = (
    @assert psm in ["C3", "C4", "C3Cytochrome"] "Photosynthesis model must be within [C3, C4, C3CytochromeModel]";

    # determine how many layers of roots
    _n_root  = 0;
    _r_index = Int[];
    for i in eachindex(z_soil)
        _z = z_soil[i];
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
    for i in eachindex(z_air)
        _z = z_air[i]
        if _z < zc
            _n_canopy += 1;
            push!(_c_index, i);
        elseif _z >= zc
            break
        end;
    end;

    # create evenly distributed root system for now
    _roots = RootHydraulics{FT}[];
    for i in _r_index
        _Δh = abs(z_soil[i+1] + z_soil[i]) / 2;
        _rt = RootHydraulics{FT}(area=1/_n_root, k_x=25/_n_root, Δh=_Δh);
        push!(_roots, _rt);
    end;

    # create leaves
    _leaves = [Leaf{FT}(psm) for i in 1:_n_canopy];
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
                zeros(FT,_n_root),  # _ks
                zeros(FT,_n_root),  # _ps
                zeros(FT,_n_root),  # _qs
    )
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-25: SPAC system for monospecies palm
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
    ROOTS::Vector{RootHydraulics{FT}}
    "Corresponding soil layer per root layer"
    ROOTS_INDEX::Vector{Int}
    "Trunk hydraulic system"
    TRUNK::StemHydraulics{FT}

    # caches to speed up calculations
    "Conductances for each root layer at given flow"
    _ks::Array{FT,1}
    "Pressure for each root layer at given flow"
    _ps::Array{FT,1}
    "Flow rate per root layer"
    _qs::Array{FT,1}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-May-25: add constructor function
#
#######################################################################################################################################################################################################
"""

    MonoPalmSPAC{FT}(psm::String; zr::Number = -1, zt::Number = 10, zc::Number = 12, z_soil::Vector = collect(0:-0.25:-2), z_air::Vector = collect(0:0.2:13)) where {FT<:AbstractFloat}

Construct a SPAC system for monospecies palm system, given
- `psm` Photosynthesis model, must be C3 or C3Cytochrome
- `zr` Maximal root depth (negative value)
- `zc` Maximal canopy height (positive value)
- `z_soil` Array of soil layer boundaries starting from 0
- `z_air` Array of air layer boundaries starting from 0

---
# Examples
```julia
spac = MonoPalmSPAC{Float64}();
spac = MonoPalmSPAC{Float64}(zr = -1, zt = 11, zc = 1, z_soil = collect(0:-0.1:-2), z_air = collect(0:0.2:13));
```
"""
MonoPalmSPAC{FT}(psm::String; zr::Number = -1, zt::Number = 10, zc::Number = 12, z_soil::Vector = collect(0:-0.25:-2), z_air::Vector = collect(0:0.2:13)) where {FT<:AbstractFloat} = (
    @assert psm in ["C3", "C3Cytochrome"] "Photosynthesis model must be within [C3, C3CytochromeModel]";

    # determine how many layers of roots
    _n_root  = 0;
    _r_index = Int[];
    for i in eachindex(z_soil)
        _z = z_soil[i];
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
    for i in 1:(length(z_air)-1)
        _z0 = z_air[i];
        _z1 = z_air[i+1];
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
    _roots = RootHydraulics{FT}[];
    for i in _r_index
        _Δh = abs(z_soil[i+1] + z_soil[i]) / 2;
        _rt = RootHydraulics{FT}(area=1/_n_root, k_x=25/_n_root, Δh=_Δh);
        push!(_roots, _rt);
    end;

    # create trunk
    _trunk = StemHydraulics{FT}(Δh=zt, Δl=zt);

    # create leaves
    _leaves = [Leaf{FT}(psm) for i in 1:_n_canopy];
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
                zeros(FT,_n_root),  # _ks
                zeros(FT,_n_root),  # _ps
                zeros(FT,_n_root),  # _qs
    )
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-25: SPAC system for monospecies tree
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
    BRANCHES::Vector{StemHydraulics{FT}}
    "Leaf hydrualic system"
    LEAVES::Vector{Leaf{FT}}
    "Corresponding air layer per canopy layer"
    LEAVES_INDEX::Vector{Int}
    "Number of canopy layers"
    N_CANOPY::Int
    "Number of root layers"
    N_ROOT::Int
    "Root hydraulic system"
    ROOTS::Vector{RootHydraulics{FT}}
    "Corresponding soil layer per root layer"
    ROOTS_INDEX::Vector{Int}
    "Trunk hydraulic system"
    TRUNK::StemHydraulics{FT}

    # caches to speed up calculations
    "Conductances for each root layer at given flow"
    _ks::Array{FT,1}
    "Pressure for each root layer at given flow"
    _ps::Array{FT,1}
    "Flow rate per root layer"
    _qs::Array{FT,1}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-May-25: add constructor function
#
#######################################################################################################################################################################################################
"""

    MonoTreeSPAC{FT}(psm::String; zr::Number = -1, zt::Number = 10, zc::Number = 12, z_soil::Vector = collect(0:-0.25:-2), z_air::Vector = collect(0:0.2:13)) where {FT<:AbstractFloat}

Construct a SPAC system for monospecies tree system, given
- `psm` Photosynthesis model, must be C3, C4, or C3Cytochrome (note: there are C4 shrubs)
- `zr` Maximal root depth (negative value)
- `zc` Maximal canopy height (positive value)
- `z_soil` Array of soil layer boundaries starting from 0
- `z_air` Array of air layer boundaries starting from 0

---
# Examples
```julia
spac = MonoTreeSPAC{Float64}();
spac = MonoTreeSPAC{Float64}(zr = -1, zt = 11, zc = 1, z_soil = collect(0:-0.1:-2), z_air = collect(0:0.2:13));
```
"""
MonoTreeSPAC{FT}(psm::String; zr::Number = -1, zt::Number = 10, zc::Number = 12, z_soil::Vector = collect(0:-0.25:-2), z_air::Vector = collect(0:0.2:13)) where {FT<:AbstractFloat} = (
    @assert psm in ["C3", "C4", "C3Cytochrome"] "Photosynthesis model must be within [C3, C4, C3CytochromeModel]";

    # determine how many layers of roots
    _n_root  = 0;
    _r_index = Int[];
    for i in eachindex(z_soil)
        _z = z_soil[i];
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
    for i in 1:(length(z_air)-1)
        _z0 = z_air[i];
        _z1 = z_air[i+1];
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
    _roots = RootHydraulics{FT}[];
    for i in _r_index
        _Δh = abs(z_soil[i+1] + z_soil[i]) / 2;
        _rt = RootHydraulics{FT}(area=1/_n_root, k_x=25/_n_root, Δh=_Δh);
        push!(_roots, _rt);
    end;

    # create trunk
    _trunk = StemHydraulics{FT}(Δh=zt, Δl=zt);

    # create branches
    _branches = StemHydraulics{FT}[];
    for i in _c_index
        _Δh = (z_air[i] + max(z_air[i+1], zt)) / 2 - zt;
        _st = StemHydraulics{FT}(area=1/_n_canopy, Δh=_Δh);
        push!(_branches, _st);
    end;

    # create leaves
    _leaves = [Leaf{FT}(psm) for i in 1:_n_canopy];
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
                zeros(FT,_n_root),  # _ks
                zeros(FT,_n_root),  # _ps
                zeros(FT,_n_root),  # _qs
    )
);
