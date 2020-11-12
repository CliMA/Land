# # Big Leaf Model

## load packages
using CanopyLayers
using PlotPlants
FT = Float32;
#------------------------------------------------------------------------------




# ## Initialization
# Besides the individual functions to initialize parameters for `CanopyLayers`,
#     a general function is provided to initialize all the parameters directly.
angles, can, can_opt, can_rad, in_rad, leaves, rt_con, rt_dim, soil, wls =
    initialize_rt_module(FT; nLayer=20, LAI3);
#------------------------------------------------------------------------------




# ## Steps
# ### 1. Update canopy optical properties (required)
canopy_geometry!(can, angles, can_opt, rt_con);
#------------------------------------------------------------------------------

# ### 2. Update scattering coefficients (required)
canopy_matrices!(leaves, can_opt);
#------------------------------------------------------------------------------

# ### 3. Simulate short wave simulation (required)
short_wave!(can, can_opt, can_rad, in_rad, soil, rt_con);
#------------------------------------------------------------------------------

# ### 4. Update integrated radiation fluxes (required for photosynthesis)
canopy_fluxes!(can, can_opt, can_rad, in_rad, soil, leaves, wls, rt_con);
#------------------------------------------------------------------------------

# ### 5. Update SIF related spectrum (required for SIF)
SIF_fluxes!(leaves, can_opt, can_rad, can, soil, wls, rt_con, rt_dim);
#------------------------------------------------------------------------------

# ### 6. Update thermo fluxes (required for leaf energy budget)
thermal_fluxes!(leaves, can_opt, can_rad, can, soil, [FT(400.0)], wls);
#------------------------------------------------------------------------------
