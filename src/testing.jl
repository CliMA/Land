# Functions and Types to be tested

# abstract typpe of plant, may contain (C3) Tree, Palm, and Grass
abstract type Plant end

# useful for tree and bamboo (with trunk and "branch")
Base.@kwdef mutable struct Tree{FT<:AbstractFloat} <: Plant
    # structre information
    "Root Layers"
    n_root  ::Int
    "Canopy Layers"
    n_canopy::Int

    # Arrays of roots and leaves
    "Roots system"
    roots ::Array{RootHydraulics{FT},1}
    "Trunk"
    trunk ::StemHydraulics{FT}
    "Branch system"
    branch::Array{StemHydraulics{FT},1}
    "Leaves"
    leaves::Array{LeafHydraulics{FT},1}

    # Root and canopy index in Soil and Atmosphere
    "Corresponding soil layer per root layer"
    root_index_in_soil ::Array{Int,1}
    "Corresponding air layer per canopy layer"
    canopy_index_in_air::Array{Int,1}
end

# useful for palm-like tree (with trunk, no "branch")
Base.@kwdef mutable struct Palm{FT<:AbstractFloat} <: Plant
    # structre information
    "Root Layers"
    n_root  ::Int
    "Canopy Layers"
    n_canopy::Int

    # Arrays of roots and leaves
    "Roots system"
    roots ::Array{RootHydraulics{FT},1}
    "Trunk"
    trunk ::StemHydraulics{FT}
    "Leaves"
    leaves::Array{LeafHydraulics{FT},1}

    # Root and canopy index in Soil and Atmosphere
    "Corresponding soil layer per root layer"
    root_index_in_soil ::Array{Int,1}
    "Corresponding air layer per canopy layer"
    canopy_index_in_air::Array{Int,1}
end

# useful for short grass (assume no stem)
Base.@kwdef mutable struct Grass{FT<:AbstractFloat} <: Plant
    # structre information
    "Root Layers"
    n_root  ::Int
    "Canopy Layers"
    n_canopy::Int

    # Arrays of roots and leaves
    "Roots system"
    roots ::Array{RootHydraulics{FT},1}
    "Leaves"
    leaves::Array{LeafHydraulics{FT},1}

    # Root and canopy index in Soil and Atmosphere
    "Corresponding soil layer per root layer"
    root_index_in_soil ::Array{Int,1}
    "Corresponding air layer per canopy layer"
    canopy_index_in_air::Array{Int,1}
end




# just create the default tree
function create_tree(
            h_root::FT,
            h_trunk::FT,
            h_canopy::FT,
            soil_bounds::Array{FT,1},
            air_bounds::Array{FT,1}
            ) where {FT<:AbstractFloat}
    # determine how many layers in roots and canopy
    _n_root  = 0;
    _r_index = Int[]
    for i in eachindex(soil_bounds)
        _z = soil_bounds[i];
        if _z > h_root
            _n_root += 1;
            push!(_r_index, i);
        else
            break
        end
    end

    _n_canopy = 0;
    _c_index  = Int[];
    for i in eachindex(air_bounds)
        _z = air_bounds[i]
        if (_z >= h_trunk) && (_z < h_canopy)
            _n_canopy += 1;
            push!(_c_index, i);
        elseif _z >= h_canopy
            break
        end
    end

    # create evenly distributed root system for now
    _roots = RootHydraulics{FT}[];
    for i in _r_index
        _Δh = (soil_bounds[i+1] + soil_bounds[i]) / 2;
        _rt = RootHydraulics{FT}(
                    area=0.1/_n_root,
                    k_max=2.5/_n_root,
                    k_rhiz=1e14/_n_root,
                    Δh=_Δh);
        push!(_roots, _rt);
    end

    # create Trunk
    _k_trunk = h_canopy / h_trunk * 5;
    _trunk   = StemHydraulics{FT}(
                    k_max=_k_trunk,
                    Δh=h_trunk);

    # create Branches
    _k_branch = h_canopy / (h_canopy - h_trunk) * 5;
    _branch   = StemHydraulics{FT}[];
    for i in _c_index
        _Δh = (air_bounds[i] + max(air_bounds[i+1], h_trunk)) / 2 - h_trunk;
        _st = StemHydraulics{FT}(
                    area=0.1/_n_canopy,
                    k_max=_k_branch/_n_canopy,
                    Δh=_Δh);
        push!(_branch, _st);
    end

    # create leaves
    _leaves = LeafHydraulics{FT}[LeafHydraulics{FT}() for i in 1:_n_canopy];

    # return plant
    return Tree{FT}(_n_root,
                    _n_canopy,
                    _roots,
                    _trunk,
                    _branch,
                    _leaves,
                    _r_index,
                    _c_index)
end




# just create a palm-like tree
function create_palm(
            h_root::FT,
            h_trunk::FT,
            h_canopy::FT,
            soil_bounds::Array{FT,1},
            air_bounds::Array{FT,1}
            ) where {FT<:AbstractFloat}
    # determine how many layers in roots and canopy
    _n_root  = 0;
    _r_index = Int[]
    for i in eachindex(soil_bounds)
        _z = soil_bounds[i];
        if _z > h_root
            _n_root += 1;
            push!(_r_index, i);
        else
            break
        end
    end

    _n_canopy = 0;
    _c_index  = Int[];
    for i in eachindex(air_bounds)
        _z = air_bounds[i]
        if (_z >= h_trunk) && (_z < h_canopy)
            _n_canopy += 1;
            push!(_c_index, i);
        elseif _z >= h_canopy
            break
        end
    end

    # create evenly distributed root system for now
    _roots = RootHydraulics{FT}[];
    for i in _r_index
        _Δh = (soil_bounds[i+1] + soil_bounds[i]) / 2;
        _rt = RootHydraulics{FT}(
                    area=0.1/_n_root,
                    k_max=2.5/_n_root,
                    k_rhiz=1e14/_n_root,
                    Δh=_Δh);
        push!(_roots, _rt);
    end

    # create Trunk
    _trunk   = StemHydraulics{FT}(
                    k_max=5,
                    Δh=h_trunk);

    # create leaves
    _leaves = LeafHydraulics{FT}[LeafHydraulics{FT}() for i in 1:_n_canopy];

    # return plant
    return Palm{FT}(_n_root,
                    _n_canopy,
                    _roots,
                    _trunk,
                    _leaves,
                    _r_index,
                    _c_index)
end
