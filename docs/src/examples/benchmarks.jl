# # Examples

## load packages
using BenchmarkTools
using WaterPhysics
FT = Float32;
#------------------------------------------------------------------------------

## define the variables
rand_r = (rand(FT) + 20) * FT(1e-6);
rand_T = rand(FT) + 298;
rand_α = rand(FT) * 50;
rand_Ψ = rand(FT) - 3;

gas_air = TraceGasAir{FT}();
gas_CO₂ = TraceGasCO₂{FT}();
#------------------------------------------------------------------------------




# ## Capillary pressure
@info "Benchmarking capillary_pressure functions...";
@btime capillary_pressure($rand_r, $rand_T);
@btime capillary_pressure($rand_r, $rand_T, $rand_α);
#------------------------------------------------------------------------------

# ## Diffusive coefficient
@info "Benchmarking diffusive_coefficient functions...";
@btime diffusive_coefficient($rand_T, $gas_CO₂, $gas_air);
@btime relative_diffusive_coefficient($rand_T);
#------------------------------------------------------------------------------

# ## Latent heat of evaporation
@info "Benchmarking latent_heat_vapor functions...";
@btime latent_heat_vapor($rand_T);
#------------------------------------------------------------------------------

# ## Saturation vapor pressure
@info "Benchmarking saturation_vapor_pressure functions...";
@btime pressure_correction($rand_T, $rand_Ψ);
@btime saturation_vapor_pressure($rand_T);
@btime saturation_vapor_pressure($rand_T, $rand_Ψ);
@btime saturation_vapor_pressure_slope($rand_T);
@btime saturation_vapor_pressure_slope($rand_T, $rand_Ψ);
#------------------------------------------------------------------------------

# ## Surface tension
@info "Benchmarking surface_tension functions...";
@btime surface_tension($rand_T);
@btime relative_surface_tension($rand_T);
#------------------------------------------------------------------------------

# ## Viscosity
@info "Benchmarking viscosity functions...";
@btime viscosity($rand_T);
@btime relative_viscosity($rand_T);
#------------------------------------------------------------------------------
