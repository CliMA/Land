###############################################################################
#
# Packages notes
# ConstrainedRootSolvers is registered, install it like normal Packages
# CanopyLayers is registered, and will be available 3 days later
# For now install it using
# pkg> add https://github.com/Yujie-W/CanopyLayers.jl
#
###############################################################################
using CanopyLayers
using ConstrainedRootSolvers
using PyPlot
using Statistics

# initialize the parameters and functions
FT   = Float64;
CAN  = create_canopy_rt(FT, nLayer=20);
WLS  = create_wave_length(FT);
DIM  = create_rt_dims(CAN, WLS);
LEAF = create_leaf_bios(FT, DIM);
fluspect!(LEAF, WLS);

MASK_PAR = (WLS.WL .> 400) .* (WLS.WL .<  700);
MASK_NIR = (WLS.WL .> 700) .* (WLS.WL .< 2500);

REF_REF_PAR = 0.0735;
REF_REF_NIR = 0.3912;
REF_TRA_PAR = 0.0566;
REF_TRA_NIR = 0.4146;

REFS = [REF_REF_PAR, REF_REF_NIR, REF_TRA_PAR, REF_TRA_NIR];
FITS = zeros(FT,4);




###############################################################################
#
# The fitting part
#
###############################################################################
# fit the values
function spectrum_diff(leaf::LeafBios{FT}, vars::Array{FT,1}, output::Array{FT,1}) where {FT<:AbstractFloat}
    leaf.N   = vars[1];
    leaf.Cab = vars[2];
    leaf.Car = vars[3];
    leaf.Ant = vars[4];
    leaf.Cs  = vars[5];
    leaf.Cw  = vars[6];
    leaf.Cm  = vars[7];
    leaf.Cx  = vars[8];
    leaf.fqe = vars[9];

    fluspect!(leaf, WLS);

    output[1] = mean(leaf.ρ_SW[MASK_PAR]);
    output[2] = mean(leaf.ρ_SW[MASK_NIR]);
    output[3] = mean(leaf.τ_SW[MASK_PAR]);
    output[4] = mean(leaf.τ_SW[MASK_NIR]);

    return -sum( (output .- REFS) .^ 2 )
end

x_inis = FT[1.5,  40, 10,  8, 0.1, 0.015, 0.012,   0,   1];
x_mins = FT[  1,   0,  0,  0,   0,     0,     0,   0,   0];
x_maxs = FT[  3, 100, 30, 40,   1,  0.05,   0.5,   1,   1];
Δ_inis = FT[0.2,  10,  3,  4, 0.1, 0.005,  0.05, 0.1, 0.1];
Δ_tole = Δ_inis .* 0.0009;

ms    = ReduceStepMethodND{FT}(x_mins = x_mins,
                               x_maxs = x_maxs,
                               x_inis = x_inis,
                               Δ_inis = Δ_inis);
st    = SolutionToleranceND{FT}(Δ_tole, 50);
_f(x) = spectrum_diff(LEAF, x, FITS);
sol   = find_peak(_f, ms, st);

@show sol;
@show spectrum_diff(LEAF, sol, FITS);




###############################################################################
#
# The plotting part
#
###############################################################################
figure(1, figsize=(10,5));
clf();
subplot(1,2,1);
plot(WLS.WL, LEAF.τ_SW, "r--", label="Leaf Tran");
plot(WLS.WL, LEAF.ρ_SW, "b--", label="Leaf Refl");
legend(loc="lower right");
title("After Fitting");

subplot(1,2,2);
plot(REFS, FITS, "ko", mfc="none", markersize=15);
plot([0,0.5], [0, 0.5], "k:");
title("After Fitting");
TEXTS = ["RP", "RN", "TP", "TN"];
for i in 1:4
    text(REFS[i], FITS[i], TEXTS[i], ha="center", va="center");
end
