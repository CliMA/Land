###############################################################################
#
# Initialize the RT module
#
###############################################################################
"""
    initialize_rt_module(FT; nLayer::Int=20, LAI::FT=FT(3.0))

Initialize the RT module so as to interface with Plant, given
- `nLayer` Number of canopy layers
- `LAI` Leaf area index

Note it here that the struct_canopy is used in CanopyRT module.
This function initializes `canopy_rt`, `canOpt_rt`, `canRad_rt`, `arrayOfLeaves` as local variables rather than global variables as in CanopyRT module.
See CanopyRT module for further operations on these variables.
"""
function initialize_rt_module(FT;
            nLayer::Int = 20,
            LAI::Number = FT(3.0)
)
    # 1. create wl_set and canopy_rt, which are rt_dims independent
    canopy_rt = create_canopy_rt(FT, nLayer=nLayer);
    wl_set    = create_wave_length(FT);

    # 2. create rt_dims from wl_set and canopy_rt
    rt_dim    = create_rt_dims(canopy_rt, wl_set);

    # 3. create canRad_rt, canOpt_rt, and etc from rt_dim and wl_set
    canRad_rt = create_canopy_rads(FT, rt_dim);
    canOpt_rt = create_canopy_opticals(FT, rt_dim);
    sunRad_rt = create_incoming_radiation(wl_set);
    soil      = create_soil_opticals(wl_set);
    angles    = SolarAngles{FT}();
    rt_con    = create_rt_container(FT, rt_dim);

    # Create an array of standard leaves
    arrayOfLeaves = [create_leaf_bios(FT, rt_dim) for i in 1:nLayer];
    for i in 1:nLayer
        fluspect!(arrayOfLeaves[i], wl_set);
    end

    # Four Different steps to compute Short-Wave RT
    #println("    compute short-wave RT...")
    canopy_geometry!(canopy_rt, angles, canOpt_rt, rt_con);
    canopy_matrices!(arrayOfLeaves, canOpt_rt);
    short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil, rt_con);
    canopy_fluxes!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil, arrayOfLeaves, wl_set, rt_con);
    SIF_fluxes!(arrayOfLeaves, canOpt_rt, canRad_rt, canopy_rt, soil, wl_set, rt_con, rt_dim);

    # # Compute Long Wave (Last term is LW incoming in W m^-2)
    #println("    compute long-wave RT...")
    thermal_fluxes!(arrayOfLeaves, canOpt_rt, canRad_rt, canopy_rt, soil, [FT(400.0)], wl_set);

    return angles, arrayOfLeaves, canopy_rt, canOpt_rt, canRad_rt, soil, sunRad_rt, wl_set, rt_con
end
