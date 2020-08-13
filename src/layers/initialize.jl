###############################################################################
#
# Initialize the RT module
#
###############################################################################
"""
    initialize_rt_module(; n_layer::Int=20, LAI::FT=FT(3.0))

Initialize the RT module so as to interface with Plant, given
- `n_layer` Number of canopy layers
- `LAI` Leaf area index

Note it here that the struct_canopy is used in CanopyRT module.
This function initializes `canopy_rt`, `canOpt_rt`, `canRad_rt`, `arrayOfLeaves` as local variables rather than global variables as in CanopyRT module.
See CanopyRT module for further operations on these variables.
"""
function initialize_rt_module(; n_layer::Int=20, LAI::FT=FT(3.0)) where {FT}
    # create canopy struct
    # to signal the difference, the variables from RT module has a postfix of _rt
    canopy_rt = Canopy4RT{FT}(nLayer=n_layer, LAI=LAI);
    wl_set    = create_wave_length(FT);
    canRad_rt = CanopyRads{FT}(nWL=wl_set.nwl, nWLf=wl_set.nWlF, nIncl=length(canopy_rt.litab), nAzi=length(canopy_rt.lazitab), nLayer=canopy_rt.nLayer);
    canOpt_rt = create_canopy_opticals(FT, wl_set.nwl, canopy_rt.nLayer, length(canopy_rt.lazitab), length(canopy_rt.litab));
    sunRad_rt = create_incoming_radiation(wl_set.swl);
    soil      = create_soil_opticals(wl_set);
    angles    = SolarAngles{FT}();
    rt_con    = create_rt_container(canopy_rt, canOpt_rt, angles, soil, wl_set);

    # Create an array of standard leaves (needs to be in Module later on:
    #println("    create leaves...")
    arrayOfLeaves = [create_leaf_bios(FT, wl_set.nwl, wl_set.nWlE, wl_set.nWlF) for i in 1:canopy_rt.nLayer];
    for i in 1:canopy_rt.nLayer
        fluspect!(arrayOfLeaves[i], wl_set);
    end

    # Four Different steps to compute Short-Wave RT
    #println("    compute short-wave RT...")
    canopy_geometry!(canopy_rt, angles, canOpt_rt, rt_con);
    canopy_matrices!(arrayOfLeaves, canOpt_rt);
    short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil, rt_con);
    canopy_fluxes!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil, arrayOfLeaves, wl_set, rt_con);

    # # Compute Long Wave (Last term is LW incoming in W m^-2)
    #println("    compute long-wave RT...")
    thermal_fluxes!(arrayOfLeaves, canOpt_rt, canRad_rt, canopy_rt, soil, [FT(400.0)], wl_set);

    return angles, arrayOfLeaves, canopy_rt, canOpt_rt, canRad_rt, soil, sunRad_rt, wl_set, rt_con
end
