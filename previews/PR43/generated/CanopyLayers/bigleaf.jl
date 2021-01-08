# load packages
using CanopyLayers
using PlotPlants
FT = Float32;

lais   = collect(FT, 0.1:0.1:5.0);
lai_sl = similar(lais);
par_sl = similar(lais); par_sh = similar(lais);
rad_sl = similar(lais); rad_sh = similar(lais);
for i in eachindex(lais)
    lai_sl[i], par_sl[i], par_sh[i], rad_sl[i], rad_sh[i] =
        big_leaf_partition(lais[i], FT(30), FT(800));
end

_fig,_axes = create_canvas("vs LAI"; ncol=3);
_ax1,_ax2,_ax3 = _axes;
_ax1.plot(lais, lai_sl, "k-");
_ax2.plot(lais, par_sl, "k-", label="Sunlit");
_ax2.plot(lais, par_sh, "k:", label="Shaded");
_ax3.plot(lais, rad_sl, "k-", label="Sunlit");
_ax3.plot(lais, rad_sh, "k:", label="Shaded");
_ax2.legend(loc="upper right");
_ax3.legend(loc="upper right");
set_xlabels!(_axes, ["Leaf area index" for i in 1:3]);
set_ylabels!(_axes, ["Sunlit fraction", "PAR (μmol m⁻² s⁻¹)", "Rabs (W m⁻²)"]);
_fig.set_tight_layout(true);
_fig

angles = collect(FT, 5:5:75);
lai_sl = similar(angles);
par_sl = similar(angles); par_sh = similar(angles);
rad_sl = similar(angles); rad_sh = similar(angles);
for i in eachindex(angles)
    lai_sl[i], par_sl[i], par_sh[i], rad_sl[i], rad_sh[i] =
        big_leaf_partition(FT(3), angles[i], FT(800));
end

_fig,_axes = create_canvas("vs Zenith angle"; ncol=3);
_ax1,_ax2,_ax3 = _axes;
_ax1.plot(angles, lai_sl, "k-");
_ax2.plot(angles, par_sl, "k-", label="Sunlit");
_ax2.plot(angles, par_sh, "k:", label="Shaded");
_ax3.plot(angles, rad_sl, "k-", label="Sunlit");
_ax3.plot(angles, rad_sh, "k:", label="Shaded");
_ax2.legend(loc="upper left");
_ax3.legend(loc="upper left");
set_xlabels!(_axes, ["Zenith angle (°)" for i in 1:3]);
set_ylabels!(_axes, ["Sunlit fraction", "PAR (μmol m⁻² s⁻¹)", "Rabs (W m⁻²)"]);
_fig.set_tight_layout(true);
_fig

rads   = collect(FT, 50:50:1000);
lai_sl = similar(rads);
par_sl = similar(rads); par_sh = similar(rads);
rad_sl = similar(rads); rad_sh = similar(rads);
for i in eachindex(rads)
    lai_sl[i], par_sl[i], par_sh[i], rad_sl[i], rad_sh[i] =
        big_leaf_partition(FT(3), FT(30), rads[i]);
end

_fig,_axes = create_canvas("vs Total radiation"; ncol=3);
_ax1,_ax2,_ax3 = _axes;
_ax1.plot(rads, lai_sl, "k-");
_ax2.plot(rads, par_sl, "k-", label="Sunlit");
_ax2.plot(rads, par_sh, "k:", label="Shaded");
_ax3.plot(rads, rad_sl, "k-", label="Sunlit");
_ax3.plot(rads, rad_sh, "k:", label="Shaded");
_ax2.legend(loc="upper left");
_ax3.legend(loc="upper left");
set_xlabels!(_axes, ["Total radiation (W m⁻²)" for i in 1:3]);
set_ylabels!(_axes, ["Sunlit fraction", "PAR (μmol m⁻² s⁻¹)", "Rabs (W m⁻²)"]);
_fig.set_tight_layout(true);
_fig

dirs   = collect(FT, 0.2:0.05:0.8);
lai_sl = similar(dirs);
par_sl = similar(dirs); par_sh = similar(dirs);
rad_sl = similar(dirs); rad_sh = similar(dirs);
for i in eachindex(dirs)
    lai_sl[i], par_sl[i], par_sh[i], rad_sl[i], rad_sh[i] =
        big_leaf_partition(FT(3), FT(30), FT(800), dirs[i]);
end

_fig,_axes = create_canvas("vs Direct light fraction"; ncol=3);
_ax1,_ax2,_ax3 = _axes;
_ax1.plot(dirs, lai_sl, "k-");
_ax2.plot(dirs, par_sl, "k-", label="Sunlit");
_ax2.plot(dirs, par_sh, "k:", label="Shaded");
_ax3.plot(dirs, rad_sl, "k-", label="Sunlit");
_ax3.plot(dirs, rad_sh, "k:", label="Shaded");
_ax2.legend(loc="upper left");
_ax3.legend(loc="upper left");
set_xlabels!(_axes, ["Direct light fraction" for i in 1:3]);
set_ylabels!(_axes, ["Sunlit fraction", "PAR (μmol m⁻² s⁻¹)", "Rabs (W m⁻²)"]);
_fig.set_tight_layout(true);
_fig

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

