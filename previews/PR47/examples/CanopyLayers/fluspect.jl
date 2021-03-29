# # Fluspect

## load packages
using Land.CanopyLayers
using PlotPlants
FT = Float32;
#------------------------------------------------------------------------------




# ## Excitation wavelength
wls    = create_wave_length(FT);
can    = create_canopy_rt(FT);
rt_dim = create_rt_dims(can, wls);
leaf   = create_leaf_bios(FT, rt_dim);
fluspect!(leaf, wls);

_fig,_axes = create_canvas("Fluspect example"; ncol=3);
_ax1,_ax2,_ax3 = _axes;
_ax1.plot(wls.WL, leaf.ρ_SW, "k-", label="Reflectance");
_ax1.plot(wls.WL, leaf.τ_SW, "k:", label="Transmittance");
_ax2.contourf(wls.WLE, wls.WLF, leaf.Mb);
_ax3.contourf(wls.WLE, wls.WLF, leaf.Mf);
_ax1.legend(loc="upper right");
_ax1.set_ylim(0,0.65);
set_xlabels!(_axes, ["Wavelength (nm)", "Excitation wavelength (nm)",
                     "Excitation wavelength (nm)"], fontsize=12);
set_ylabels!(_axes, ["Reflectance or Transmittance", "SIF wavelength (nm)",
                     "SIF wavelength (nm)"], fontsize=12);
set_titles!(_axes; labels=["Spectrum", "Backward matrix", "Forward matrix"],
                   usetex=false);
_fig.set_tight_layout(true);
_fig
#------------------------------------------------------------------------------




# ## Change leaf properties
## here we change all the properties at the same time as an example
leaf.N   = 2.0;
leaf.Cab = 50.0;
leaf.Car = 15.0;
leaf.Ant = 0.1;
leaf.Cs  = 0.1;
leaf.Cw  = 0.02;
leaf.Cm  = 0.02;
leaf.Cx  = 0.1;
leaf.fqe = 0.8;
fluspect!(leaf, wls);

_fig,_axes = create_canvas("Change leaf properties"; ncol=3);
_ax1,_ax2,_ax3 = _axes;
_ax1.plot(wls.WL, leaf.ρ_SW, "k-", label="Reflectance");
_ax1.plot(wls.WL, leaf.τ_SW, "k:", label="Transmittance");
_ax2.contourf(wls.WLE, wls.WLF, leaf.Mb);
_ax3.contourf(wls.WLE, wls.WLF, leaf.Mf);
_ax1.legend(loc="upper right");
_ax1.set_ylim(0,0.65);
set_xlabels!(_axes, ["Wavelength (nm)", "Excitation wavelength (nm)",
                     "Excitation wavelength (nm)"], fontsize=12);
set_ylabels!(_axes, ["Reflectance or Transmittance", "SIF wavelength (nm)",
                     "SIF wavelength (nm)"], fontsize=12);
set_titles!(_axes; labels=["Spectrum", "Backward matrix", "Forward matrix"],
                   usetex=false);
_fig.set_tight_layout(true);
_fig
#------------------------------------------------------------------------------
