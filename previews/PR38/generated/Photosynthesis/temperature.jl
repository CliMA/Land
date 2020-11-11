# load packages
using Land.Photosynthesis
using PlotPlants
FT = Float32;

# KcTDCLM is a ArrheniusTD type struct, use it as an example here
_td = KcTDCLM(FT);
_ts = collect(FT, 273:1:323);
_ks = temperature_correction.([_td], _ts);

_fig,_axes = create_canvas("Arrhenius correction without deactivation");
_ax1 = _axes[1];
_ax1.plot(_ts .- 273.15, _ks, "k-");
_ax1.set_xlabel("Leaf temperature (°C)");
_ax1.set_ylabel("Relative to 25 °C");
_fig

# VcmaxTDCLM is a ArrheniusPeakTD type struct, use it as an example here
_td = VcmaxTDCLM(FT);
_ts = collect(FT, 273:1:323);
_ks = temperature_correction.([_td], _ts);

_fig,_axes = create_canvas("Arrhenius correction with deactivation");
_ax1 = _axes[1];
_ax1.plot(_ts .- 273.15, _ks, "k-");
_ax1.set_xlabel("Leaf temperature (°C)");
_ax1.set_ylabel("Relative to 25 °C");
_fig

_td = Q10TD{FT}(1.0, 273.15, 1.7);
_ts = collect(FT, 273:1:323);
_ks = temperature_correction.([_td], _ts);

_fig,_axes = create_canvas("Q10 correction");
_ax1 = _axes[1];
_ax1.plot(_ts .- 273.15, _ks, "k-");
_ax1.set_xlabel("Leaf temperature (°C)");
_ax1.set_ylabel("Relative to 25 °C");
_fig

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

