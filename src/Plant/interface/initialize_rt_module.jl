"""
    initialize_rt_module(n_layer)

Initialize the RT module so as to interface with Plant, given
- `n_layer` Number of canopy layers
- `LAI` Leaf area index

Note it here that the struct_canopy is used in CanopyRT module.
This function initializes `canopy_rt`, `canOpt_rt`, `canRad_rt`, `arrayOfLeaves` as local variables rather than global variables as in CanopyRT module.
See CanopyRT module for further operations on these variables.
"""
function initialize_rt_module(; n_layer::Int=10, LAI::FT=FT(3.0)) where {FT}
    # create canopy struct
    # to signal the difference, the variables from RT module has a postfix of _rt
    canopy_rt = struct_canopy{FT, n_layer, LAI}()
    canRad_rt = struct_canopyRadiation{FT, CanopyRT.nwl, CanopyRT.nWlF, length(CanopyRT.litab), length(canopy_rt.lazitab), canopy_rt.nlayers}()
    canOpt_rt = create_canopyOpt(FType=FT,nWL=CanopyRT.nwl,nLayers=canopy_rt.nlayers, nAzi=length(canopy_rt.lazitab), nIncl=length(CanopyRT.litab))

    # Create an array of standard leaves (needs to be in Module later on:
    println("    create leaves...")
    arrayOfLeaves = Array{leafbio{FT, CanopyRT.nwl, length(CanopyRT.wle), length(CanopyRT.wlf), length(CanopyRT.wle)*length(CanopyRT.wlf)}, 1}(undef, canopy_rt.nlayers)
    for i = 1:canopy_rt.nlayers
        arrayOfLeaves[i] = leafbio{FT, CanopyRT.nwl, length(CanopyRT.wle), length(CanopyRT.wlf),length(CanopyRT.wle)*length(CanopyRT.wlf)}()
        fluspect!(arrayOfLeaves[i], CanopyRT.optis)
    end

    # Four Different steps to compute Short-Wave RT
    println("    compute short-wave RT...")
    computeCanopyGeomProps!(canopy_rt, CanopyRT.angles, canOpt_rt)
    computeCanopyMatrices!(arrayOfLeaves, canOpt_rt);
    RTM_SW!(canopy_rt, canOpt_rt, canRad_rt, CanopyRT.sunRad, CanopyRT.soil);
    deriveCanopyFluxes!(canopy_rt, canOpt_rt, canRad_rt, CanopyRT.sunRad, CanopyRT.soil, arrayOfLeaves);

    # # Compute Long Wave (Last term is LW incoming in W m^-2)
    println("    compute long-wave RT...")
    computeThermalFluxes!(arrayOfLeaves, canOpt_rt, canRad_rt, canopy_rt, CanopyRT.soil, [FT(400.0)]);

    return canopy_rt, canOpt_rt, canRad_rt, arrayOfLeaves
end