# # Examples

## load packages
using PyPlot
using WaterPhysics
FT = Float32;

## define preview_data function
function preview_data(
            xs::Array{FT,1},
            ys::Array{FT,1},
            title::String,
            xlab::String,
            ylab::String
) where {FT<:AbstractFloat}
    fig = figure(title, figsize=(4,3));
    ax0 = fig.add_subplot(1,1,1);
    ax0.plot(xs, ys);
    ax0.set_xlabel(xlab);
    ax0.set_ylabel(ylab);
    fig.set_tight_layout(true);

    return fig
end
#------------------------------------------------------------------------------




# ## Capilary pressure
# ### As a function of radius
Rs = collect(FT, 1:0.5:30) .* FT(1e-6);
Pc = capillary_pressure.(Rs, FT(298.15));

preview_data(Rs .* FT(1e6), Pc .* FT(1e-3),
             "Pc ~ Radius", "Radius (μm)", "Capillary pressure (kPa)")
#------------------------------------------------------------------------------

# ### As a function of contact angle
αs = collect(FT, 0:1:90);
Pc = capillary_pressure.(FT(10 * 1e-6), FT(298.15), αs);

preview_data(αs, Pc .* FT(1e-3),
             "Pc ~ α", "Contact angle (°)", "Capillary pressure (kPa)")
#------------------------------------------------------------------------------

# ### As a function of temperature (due to surface tension)
Ts = collect(FT, 0:1:50) .+ FT(273.15);
Pc = capillary_pressure.(FT(10 * 1e-6), Ts);

preview_data(Ts .- FT(273.15), Pc .* FT(1e-3),
             "Pc ~ T", "Temperature (°C)", "Capillary pressure (kPa)")
#------------------------------------------------------------------------------




# ## Diffusive coefficient
# ### Relative coefficient as a function of temperature
Ts = collect(FT, 0:1:50) .+ FT(273.15);
Dr = relative_diffusive_coefficient.(Ts);

preview_data(Ts .- FT(273.15), Dr,
             "rDc ~ T", "Temperature (°C)", "Relative diffusive coefficient")
#------------------------------------------------------------------------------

# ### Relative coefficient as a function of temperature
Ts = collect(FT, 0:1:50) .+ FT(273.15);
Ds = diffusive_coefficient.(Ts, [TraceGasCO₂()], [TraceGasAir()]);

preview_data(Ts .- FT(273.15), Ds,
             "Dc ~ T", "Temperature (°C)", "Diffusive coefficient (m² s⁻¹)")
#------------------------------------------------------------------------------




# ## Latent heat
# ### λ of liquid water as a function of temperature
Ts = collect(FT, 0:1:50) .+ FT(273.15);
λs = latent_heat_vapor.(Ts);

preview_data(Ts .- FT(273.15), λs,
             "λ ~ T", "Temperature (°C)", "Latent heat (J kg⁻¹)")
#------------------------------------------------------------------------------




# ## Saturation vapor pressure
# ### Saturated water vapor pressure at flat surface
Ts = collect(FT, 0:1:50) .+ FT(273.15);
Ps = saturation_vapor_pressure.(Ts);

preview_data(Ts .- FT(273.15), Ps,
             "Ps ~ T", "Temperature (°C)", "Saturation vapor pressure (Pa)")
#------------------------------------------------------------------------------

# ### Saturated water vapor pressure at curved surface (Kelvin equation)
Ψs = collect(FT, -5:0.1:5) .* FT(1e6);
Ps = saturation_vapor_pressure.(FT(298.15), Ψs);

preview_data(Ψs .* FT(1e-6), Ps,
             "Ps ~ Ψs",
             "Water potential (MPa)", "Saturation vapor pressure (Pa)")
#------------------------------------------------------------------------------




# ## Surface tension
# ### Relative γ of water-air interface as a function of temperature
Ts = collect(FT, 0:1:50) .+ FT(273.15);
γr = relative_surface_tension.(Ts);

preview_data(Ts .- FT(273.15), γr,
             "rγ ~ T", "Temperature (°C)", "Relative surface tension")
#------------------------------------------------------------------------------

# ### γ of water-air interface as a function of temperature
Ts = collect(FT, 0:1:50) .+ FT(273.15);
γs = surface_tension.(Ts);

preview_data(Ts .- FT(273.15), γs,
             "γ ~ T", "Temperature (°C)", "Surface tension (N m⁻¹)")
#------------------------------------------------------------------------------




# ## Viscosity
# ### Relative μ of water as a function of temperature
Ts = collect(FT, 0:1:50) .+ FT(273.15);
μr = relative_viscosity.(Ts);

preview_data(Ts .- FT(273.15), μr,
             "rμ ~ T", "Temperature (°C)", "Relative viscosity")
#------------------------------------------------------------------------------

# ### Relative μ of water as a function of temperature
Ts = collect(FT, 0:1:50) .+ FT(273.15);
μs = viscosity.(Ts);

preview_data(Ts .- FT(273.15), μs,
             "μ ~ T", "Temperature (°C)", "Viscosity (Pa s)")
#------------------------------------------------------------------------------
