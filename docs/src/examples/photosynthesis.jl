# # Temperature dependencies

## load packages
using Photosynthesis
using PlotPlants
FT = Float32;
#------------------------------------------------------------------------------




# ## Simple example
## define photosynthesis system and leaf (C3 and C4), and envir
c3_set = C3CLM(FT);
c4_set = C4CLM(FT);
leaf_3 = Leaf{FT}(APAR=1000);
leaf_4 = Leaf{FT}(APAR=1000);
envir  = AirLayer{FT}();

## define leaf temperature, total leaf conductance to CO₂, and internal CO₂
T   = FT(300);
glc = FT(0.1);
p_i = rand(FT) + 20;

## remember to update the temperature dependencies when temperature changes
println("initialize temperature dependencies");
leaf_temperature_dependence!(c3_set, leaf_3, envir, T);
leaf_temperature_dependence!(c4_set, leaf_4, envir, T);

println("calculate photosynthesis from known internal CO₂ partial pressure");
leaf_photo_from_pi!(c3_set, leaf_3, p_i);
leaf_photo_from_pi!(c4_set, leaf_4, p_i);
@show leaf_3.An;
@show leaf_4.An;

println("calculate photosynthesis from known leaf conductance to CO₂");
leaf_photo_from_glc!(c3_set, leaf_3, envir, glc);
leaf_photo_from_glc!(c4_set, leaf_4, envir, glc);
@show leaf_3.An;
@show leaf_4.An;
#------------------------------------------------------------------------------




# ## A-Ci curve
# Here we show an example of the A-Ci curves for C3 and C4 photosynthesis. As
#     stomatal conductance may differ when other environmental conditions
#     change, we leave the examples of photosynthesis responses to the
#     environment for `StomataModels` package.
## temperature not changing, no temperature dependencies update required
_p3 = collect(FT, 5:1:200);
_p4 = collect(FT, 0:0.1:15.01);
_a3 = similar(_p3);
_a4 = similar(_p4);
for i in eachindex(_p3)
    leaf_photo_from_pi!(c3_set, leaf_3, _p3[i]); _a3[i] = leaf_3.An;
end
for i in eachindex(_p4)
    leaf_photo_from_pi!(c4_set, leaf_4, _p4[i]); _a4[i] = leaf_4.An;
end

_fig,_axes = create_canvas("A-Ci curve"; ncol=2);
_ax1,_ax2 = _axes;
_ax1.plot(_p3, _a3, "k-", label="C3");
_ax2.plot(_p4, _a4, "k-", label="C4");
_ax1.set_xlabel("Leaf internal CO₂ (Pa)");
_ax1.set_ylabel("Anet (μmol m⁻² s⁻¹)");
_ax2.set_xlabel("Leaf internal CO₂ (Pa)");
_ax1.legend(loc="lower right");
_ax2.legend(loc="lower right");
_fig
#------------------------------------------------------------------------------
