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

    MonoGrassSPAC{FT}(psm::String; z_root::Number = -0.2, z_canopy::Number = 0.5, soil_bounds::Vector = collect(0:-0.1:-1), air_bounds::Vector = collect(0:0.05:1)) where {FT<:AbstractFloat}

Construct a SPAC system for monospecies grass system, given
- `psm` Photosynthesis model, must be C3, C4 cor C3Cytochrome
- `z_root` Maximal root depth (negative value)
- `z_canopy` Maximal canopy height (positive value)
- `soil_bounds` Array of soil layer boundaries starting from 0
- `air_bounds` Array of air layer boundaries starting from 0

---
# Examples
```julia
spac = MonoGrassSPAC{Float64}();
spac = MonoGrassSPAC{Float64}(z_root = -0.3, z_canopy = 1, soil_bounds = collect(0:-0.1:-1), air_bounds = collect(0:0.05:1.01));
```
"""
MonoGrassSPAC{FT}(psm::String; z_root::Number = -0.2, z_canopy::Number = 0.5, soil_bounds::Vector = collect(0:-0.1:-1), air_bounds::Vector = collect(0:0.05:1)) where {FT<:AbstractFloat} = (
    @assert psm in ["C3", "C4", "C3Cytochrome"] "Photosynthesis model must be within [C3, C4, C3CytochromeModel]";

    # determine how many layers of roots
    _n_root  = 0;
    _r_index = Int[];
    for i in eachindex(soil_bounds)
        _z = soil_bounds[i];
        if _z > z_root
            _n_root += 1;
            push!(_r_index, i);
        else
            break
        end;
    end;

    # determine how many layers of canopy
    _n_canopy = 0;
    _c_index  = Int[];
    for i in eachindex(air_bounds)
        _z = air_bounds[i]
        if _z < z_canopy
            _n_canopy += 1;
            push!(_c_index, i);
        elseif _z >= z_canopy
            break
        end;
    end;

    # create evenly distributed root system for now
    _roots = RootHydraulics{FT}[];
    for i in _r_index
        _Δh = abs(soil_bounds[i+1] + soil_bounds[i]) / 2;
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
