###############################################################################
#
# Creating Plant Hydraulic system
#
###############################################################################
"""
    create_grass(
                z_root::FT,
                z_canopy::FT,
                soil_bounds::Vector{FT},
                air_bounds::Vector{FT}
    ) where {FT<:AbstractFloat}

Create a [`GrassLikeOrganism`](@ref), given
- `z_root` Maximal root depth (negative value)
- `z_canopy` Maximal canopy height (positive value)
- `soil_bounds` Array of soil layer boundaries starting from 0
- `air_bounds` Array of air layer boundaries starting from 0
"""
function create_grass(
            z_root::FT,
            z_canopy::FT,
            soil_bounds::Vector{FT},
            air_bounds::Vector{FT}
) where {FT<:AbstractFloat}
    # determine how many layers in roots and canopy
    _n_root  = 0;
    _r_index = Int[]
    for i in eachindex(soil_bounds)
        _z = soil_bounds[i];
        if _z > z_root
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
        if _z < z_canopy
            _n_canopy += 1;
            push!(_c_index, i);
        elseif _z >= z_canopy
            break
        end
    end

    # create evenly distributed root system for now
    _roots = RootHydraulics{FT}[];
    for i in _r_index
        _Δh = abs(soil_bounds[i+1] + soil_bounds[i]) / 2;
        _rt = RootHydraulics{FT}(
                    area=1/_n_root,
                    k_max=25/_n_root,
                    k_rhiz=5e14/_n_root,
                    Δh=_Δh);
        push!(_roots, _rt);
    end

    # create leaves
    _leaves = [LeafHydraulics{FT}(area=1500/_n_canopy) for i in 1:_n_canopy];

    # return plant
    return GrassLikeOrganism{FT}(_n_root,
                                 _n_canopy,
                                 _roots,
                                 _leaves,
                                 _r_index,
                                 _c_index,
                                 zeros(FT,_n_root),
                                 zeros(FT,_n_root),
                                 zeros(FT,_n_root))
end




"""
    create_palm(z_root::FT,
                z_trunk::FT,
                z_canopy::FT,
                soil_bounds::Vector{FT},
                air_bounds::Vector{FT}
    ) where {FT<:AbstractFloat}

Create a [`PalmLikeOrganism`](@ref), given
- `z_root` Maximal root depth (negative value)
- `z_trunk` Maximal trunk height (positive value)
- `z_canopy` Maximal canopy height (positive value)
- `soil_bounds` Array of soil layer boundaries starting from 0
- `air_bounds` Array of air layer boundaries starting from 0
"""
function create_palm(
            z_root::FT,
            z_trunk::FT,
            z_canopy::FT,
            soil_bounds::Vector{FT},
            air_bounds::Vector{FT}
) where {FT<:AbstractFloat}
    # determine how many layers in roots and canopy
    _n_root  = 0;
    _r_index = Int[]
    for i in eachindex(soil_bounds)
        _z = soil_bounds[i];
        if _z > z_root
            _n_root += 1;
            push!(_r_index, i);
        else
            break
        end
    end

    _n_canopy = 0;
    _c_index  = Int[];
    for i in 1:(length(air_bounds)-1)
        _z0 = air_bounds[i];
        _z1 = air_bounds[i+1];
        if _z0 <= z_trunk < z_canopy <= _z1
            _n_canopy += 1;
            push!(_c_index, i);
            break
        elseif _z0 < z_trunk < _z1 < z_canopy
            _n_canopy += 1;
            push!(_c_index, i);
        elseif z_trunk <= _z0 < _z1 < z_canopy
            _n_canopy += 1;
            push!(_c_index, i);
        elseif z_trunk <= _z0 <= z_canopy < _z1
            _n_canopy += 1;
            push!(_c_index, i-1);
            break
        elseif z_canopy <= _z0 < _z1
            break
        end
    end

    # create evenly distributed root system for now
    _roots = RootHydraulics{FT}[];
    for i in _r_index
        _Δh = abs(soil_bounds[i+1] + soil_bounds[i]) / 2;
        _rt = RootHydraulics{FT}(
                    area=1/_n_root,
                    k_max=25/_n_root,
                    k_rhiz=5e14/_n_root,
                    Δh=_Δh);
        push!(_roots, _rt);
    end

    # create Trunk
    _trunk   = StemHydraulics{FT}(
                    k_max=50,
                    Δh=z_trunk);

    # create leaves
    _leaves = [LeafHydraulics{FT}(area=1500/_n_canopy) for i in 1:_n_canopy];

    # return plant
    return PalmLikeOrganism{FT}(_n_root,
                                _n_canopy,
                                _roots,
                                _trunk,
                                _leaves,
                                _r_index,
                                _c_index,
                                zeros(FT,_n_root),
                                zeros(FT,_n_root),
                                zeros(FT,_n_root))
end




"""
    create_tree(z_root::FT,
                z_trunk::FT,
                z_canopy::FT,
                soil_bounds::Vector{FT},
                air_bounds::Vector{FT}
    ) where {FT<:AbstractFloat}

Create a [`TreeLikeOrganism`](@ref), given
- `z_root` Maximal root depth (negative value)
- `z_trunk` Maximal trunk height (positive value)
- `z_canopy` Maximal canopy height (positive value)
- `soil_bounds` Array of soil layer boundaries starting from 0
- `air_bounds` Array of air layer boundaries starting from 0
"""
function create_tree(
            z_root::FT,
            z_trunk::FT,
            z_canopy::FT,
            soil_bounds::Vector{FT},
            air_bounds::Vector{FT}
) where {FT<:AbstractFloat}
    # determine how many layers in roots and canopy
    _n_root  = 0;
    _r_index = Int[]
    for i in eachindex(soil_bounds)
        _z = soil_bounds[i];
        if _z > z_root
            _n_root += 1;
            push!(_r_index, i);
        else
            break
        end
    end

    _n_canopy = 0;
    _c_index  = Int[];
    for i in eachindex(air_bounds)
        _z0 = air_bounds[i];
        _z1 = air_bounds[i+1];
        if _z0 <= z_trunk < z_canopy <= _z1
            _n_canopy += 1;
            push!(_c_index, i);
            break
        elseif _z0 < z_trunk < _z1 < z_canopy
            _n_canopy += 1;
            push!(_c_index, i);
        elseif z_trunk <= _z0 < _z1 < z_canopy
            _n_canopy += 1;
            push!(_c_index, i);
        elseif z_trunk <= _z0 <= z_canopy < _z1
            _n_canopy += 1;
            push!(_c_index, i-1);
            break
        elseif z_canopy <= _z0 < _z1
            break
        end
    end

    # create evenly distributed root system for now
    _roots = RootHydraulics{FT}[];
    for i in _r_index
        _Δh = abs(soil_bounds[i+1] + soil_bounds[i]) / 2;
        _rt = RootHydraulics{FT}(
                    area=1/_n_root,
                    k_max=25/_n_root,
                    k_rhiz=5e14/_n_root,
                    Δh=_Δh);
        push!(_roots, _rt);
    end

    # create Trunk
    _k_trunk = z_canopy / z_trunk * 50;
    _trunk   = StemHydraulics{FT}(
                    k_max=_k_trunk,
                    Δh=z_trunk);

    # create Branches
    _k_branch = z_canopy / (z_canopy - z_trunk) * 50;
    _branch   = StemHydraulics{FT}[];
    for i in _c_index
        _Δh = (air_bounds[i] + max(air_bounds[i+1], z_trunk)) / 2 - z_trunk;
        _st = StemHydraulics{FT}(
                    area=1/_n_canopy,
                    k_max=_k_branch/_n_canopy,
                    Δh=_Δh);
        push!(_branch, _st);
    end

    # create leaves
    _leaves = [LeafHydraulics{FT}(area=1500/_n_canopy) for i in 1:_n_canopy];

    # return plant
    return TreeLikeOrganism{FT}(_n_root,
                                _n_canopy,
                                _roots,
                                _trunk,
                                _branch,
                                _leaves,
                                _r_index,
                                _c_index,
                                zeros(FT,_n_root),
                                zeros(FT,_n_root),
                                zeros(FT,_n_root))
end
