###############################################################################
#
# Initialize the RT module
#
###############################################################################
"""
    initialize_rt_module(FT; nLayer::Int = 20, LAI::Number = FT(3))

Initialize the RT module and return the sturctures, given
- `FT` Floating number type
- `nLayer` Number of canopy layers
- `LAI` Leaf area index

This function initializes and returns
- `angles` [`SolarAngles`](@ref)
- `can` [`Canopy4RT`](@ref)
- `can_opt` [`CanopyOpticals`](@ref)
- `can_rad` [`CanopyRads`](@ref)
- `in_rad` [`IncomingRadiation`](@ref)
- `leaves` Array{[`LeafBios`](@ref),1}
- `rt_con` [`RTCache`](@ref)
- `rt_dim` [`RTDimensions`](@ref)
- `soil` [`SoilOpticals`](@ref)
- `wls` [`WaveLengths`](@ref)
"""
function initialize_rt_module(FT; nLayer::Int = 20, LAI::Number = FT(3))
    # 1. create wls and can, which are rt_dims independent
    can = create_canopy_rt(FT, nLayer=nLayer, LAI=LAI);
    wls = create_wave_length(FT);

    # 2. create rt_dims from wls and can
    rt_dim = create_rt_dims(can, wls);

    # 3. create can_rad, can_opt, and etc from rt_dim and wls
    can_rad = create_canopy_rads(FT, rt_dim);
    can_opt = create_canopy_opticals(FT, rt_dim);
    in_rad  = create_incoming_radiation(wls);
    soil    = SoilOpticals{FT}(wls);
    angles  = SolarAngles{FT}();
    rt_con  = create_rt_cache(FT, rt_dim);

    # Create an array of standard leaves
    leaves = [create_leaf_bios(FT, rt_dim) for i in 1:nLayer];
    for i in 1:nLayer
        fluspect!(leaves[i], wls);
    end

    # Four Different steps to compute Short-Wave RT
    canopy_geometry!(can, angles, can_opt, rt_con);
    canopy_matrices!(leaves, can_opt);
    short_wave!(can, can_opt, can_rad, in_rad, soil, rt_con);
    canopy_fluxes!(can, can_opt, can_rad, in_rad, soil, leaves, wls, rt_con);
    SIF_fluxes!(leaves, can_opt, can_rad, can, soil, wls, rt_con, rt_dim);

    # # Compute Long Wave (Last term is LW incoming in W m^-2)
    thermal_fluxes!(leaves, can_opt, can_rad, can, soil, [FT(400.0)], wls);

    return angles, can, can_opt, can_rad, in_rad, leaves, rt_con, rt_dim, soil, wls
end
