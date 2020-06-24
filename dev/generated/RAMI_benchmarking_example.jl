#### Use Julia Plots package and switch to plotly js option:
using PyPlot

##----------------------------------------------------------------------------

# First, we include Revise (good for debugging) and Parameters (tools for structures)

##using Revise
using Parameters
##----------------------------------------------------------------------------

# Now include the Land modules

using Land
using Land.CanopyRT
##----------------------------------------------------------------------------

##Defining all reference values for the Sparse case

const FT = Float32

wl_set = create_wl_para_set(FT)
leaf = create_leaf_bio(FT, wl_set.nwl, wl_set.nWlE, wl_set.nWlF);
canopy_rt = Canopy4RT{FT, 20, 3.0}()
canRad_rt = CanopyRadiation{FT, wl_set.nwl, wl_set.nWlF, length(canopy_rt.litab), length(canopy_rt.lazitab), canopy_rt.nlayers}()
canOpt_rt = create_canopy_optical(FT, wl_set.nwl, canopy_rt.nlayers, length(canopy_rt.lazitab), length(canopy_rt.litab); using_marray=false)
sunRad_rt = create_incoming_radiation(FT, wl_set.swl);

soil = SoilOpti{FT}(wl_set.wl, FT(0.2)*ones(FT, length(wl_set.wl)), FT[0.1], FT(290.0))
angles = SolarAngles{FT}()

arrayOfLeaves = [create_leaf_bio(FT, wl_set.nwl, wl_set.nWlE, wl_set.nWlF) for i in 1:canopy_rt.nlayers]
for i in 1:canopy_rt.nlayers
    fluspect!(arrayOfLeaves[i], wl_set)
end

RAMI_SZA = [27.,60.,83.]

RAMI_fabsRed_050_BLK =  [0.09380509999999999, 0.16259713, 0.53931207]
RAMI_frefRed_050_BLK =  [0.00330673, 0.00517598, 0.01626682]
RAMI_ftranRed_050_BLK =  [0.90288817, 0.83222689, 0.44442110999999995]

RAMI_fabsRed_050_MED =  [0.10897124, 0.17760124000000002, 0.54764719]
RAMI_frefRed_050_MED =  [0.09759354, 0.09107608, 0.06177913]
RAMI_ftranRed_050_MED =  [0.90337609, 0.83265704, 0.44469279]

RAMI_fabsRed_050_SNW =  [0.21471034, 0.28200132, 0.60564705]
RAMI_frefRed_050_SNW =  [0.7526521700000001, 0.6879087300000001, 0.37825442000000004]
RAMI_ftranRed_050_SNW =  [0.90659694, 0.83583194, 0.44718138]

function RAMI_case(LAI, soil_albedo, clumping_index)

  soil.albedo_SW[:] .=soil_albedo;
  ##Clumping index
  canopy_rt.Ω = 1.0
  ##Viewing Zenith Angle in degrees
  angles.tto=0.0
  ##Leaf Area index
  canopy_rt.LAI=LAI

  reflRed_SZA = []
  absRed_SZA = []
  transRed_SZA = []

  for SZA=0.0:1:85
    angles.tts=SZA

    fluspect!(leaf, wl_set);
    compute_canopy_geometry!(canopy_rt, angles, canOpt_rt)
    compute_canopy_matrices!(arrayOfLeaves, canOpt_rt);

   # leaf reflectance RED
   leaf.ρ_SW[28] = 0.0735
   # leaf transmittance
   leaf.τ_SW[28]= 0.0566

    ##Setting all diffuse to zero
    sunRad_rt.E_diffuse[28] = 0.0

    simulate_short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil);
    push!(reflRed_SZA, canRad_rt.alb_direct[28])
    push!(absRed_SZA, (sum(canRad_rt.netSW_shade,dims=2)[28,1].+sum(canRad_rt.netSW_sunlit,dims=2)[28,1])./(sunRad_rt.E_diffuse[28].+sunRad_rt.E_direct[28]))
    push!(transRed_SZA,  (canOpt_rt.Es_[28,end] .+ canRad_rt.E_down[28,end])./(sunRad_rt.E_diffuse[28].+sunRad_rt.E_direct[28]))
  end

  ######## Clumped Case

  reflRed_clump_SZA = []
  absRed_clump_SZA = []
  transRed_clump_SZA = []


  ##Clumping index
  canopy_rt.Ω = clumping_index

  for SZA=0.0:1:85
    angles.tts=SZA

    fluspect!(leaf, wl_set);
    compute_canopy_geometry!(canopy_rt, angles, canOpt_rt)
    compute_canopy_matrices!(arrayOfLeaves, canOpt_rt);

    simulate_short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil);
    push!(reflRed_clump_SZA, canRad_rt.alb_direct[28])
    push!(absRed_clump_SZA, (sum(canRad_rt.netSW_shade,dims=2)[28,1].+sum(canRad_rt.netSW_sunlit,dims=2)[28,1])./(sunRad_rt.E_diffuse[28].+sunRad_rt.E_direct[28]))
    push!(transRed_clump_SZA,  (canOpt_rt.Es_[28,end] .+ canRad_rt.E_down[28,end])./(sunRad_rt.E_diffuse[28].+sunRad_rt.E_direct[28]))

  end

  return reflRed_SZA,absRed_SZA,transRed_SZA,reflRed_clump_SZA,absRed_clump_SZA,transRed_clump_SZA

end;

#TODO nest those into a loop!

##Sparse case with black soil
reflRed_SZA,absRed_SZA,transRed_SZA,reflRed_clump_SZA,absRed_clump_SZA,transRed_clump_SZA=RAMI_case(0.50265, 0.0, 0.365864235);

SZA=0:1:85

figure(figsize=(10,5))

subplot(1,2,1)
plot(SZA, reflRed_SZA,label=["reflectance"])
plot(SZA, absRed_SZA ,label=["absorptance"])
plot(SZA, transRed_SZA,label=["transmittance"])
scatter(RAMI_SZA, RAMI_frefRed_050_BLK)
scatter(RAMI_SZA, RAMI_fabsRed_050_BLK)
scatter(RAMI_SZA, RAMI_ftranRed_050_BLK)
title("050_BLK - Default")
ylabel("Radiation Partitioning")
xlabel("Sun Zenith Angle")
xlim(0.0, 90.)
ylim(-0.05, 1.0)
xticks(0:20:91.)
yticks(0:0.2:1.0)

subplot(1,2,2)
plot(SZA, reflRed_clump_SZA,label=["reflectance"])
plot(SZA, absRed_clump_SZA ,label=["absorptance"])
plot(SZA,transRed_clump_SZA,label=["transmittance"])
scatter(RAMI_SZA, RAMI_frefRed_050_BLK,label="RAMI reflectance")
scatter(RAMI_SZA, RAMI_fabsRed_050_BLK,label="RAMI absorptance")
scatter(RAMI_SZA, RAMI_ftranRed_050_BLK,label="RAMI transmittance")
title("Clumping")
ylabel("Radiation Partitioning")
xlabel("Sun Zenith Angle")
xlim(0.0, 90.)
ylim(-0.05, 1.0)
xticks(0:20:91.)
yticks(0:0.2:1.0)

legend()
gcf()

##Sparse case with medium soil
reflRed_SZA,absRed_SZA,transRed_SZA,reflRed_clump_SZA,absRed_clump_SZA,transRed_clump_SZA=RAMI_case(0.50265, 0.1217, 0.365864235);

SZA=0:1:85

figure(figsize=(10,5))

subplot(1,2,1)
plot(SZA, reflRed_SZA,label=["reflectance"])
plot(SZA, absRed_SZA ,label=["absorptance"])
plot(SZA, transRed_SZA,label=["transmittance"])
scatter(RAMI_SZA, RAMI_frefRed_050_MED)
scatter(RAMI_SZA, RAMI_fabsRed_050_MED)
scatter(RAMI_SZA, RAMI_ftranRed_050_MED)
title("050_MED - Default")
ylabel("Radiation Partitioning")
xlabel("Sun Zenith Angle")
xlim(0.0, 90.)
ylim(-0.05, 1.0)
xticks(0:20:91.)
yticks(0:0.2:1.0)

subplot(1,2,2)
plot(SZA, reflRed_clump_SZA,label=["reflectance"])
plot(SZA, absRed_clump_SZA ,label=["absorptance"])
plot(SZA,transRed_clump_SZA,label=["transmittance"])
scatter(RAMI_SZA, RAMI_frefRed_050_MED,label="RAMI reflectance")
scatter(RAMI_SZA, RAMI_fabsRed_050_MED,label="RAMI absorptance")
scatter(RAMI_SZA, RAMI_ftranRed_050_MED,label="RAMI transmittance")
title("Clumping")
ylabel("Radiation Partitioning")
xlabel("Sun Zenith Angle")
xlim(0.0, 90.)
ylim(-0.05, 1.0)
xticks(0:20:91.)
yticks(0:0.2:1.0)

legend()
gcf()

##Sparse case with snowy soil
reflRed_SZA,absRed_SZA,transRed_SZA,reflRed_clump_SZA,absRed_clump_SZA,transRed_clump_SZA=RAMI_case(0.50265, 0.9640, 0.365864235);

SZA=0:1:85

figure(figsize=(10,5))

subplot(1,2,1)
plot(SZA, reflRed_SZA,label=["reflectance"])
plot(SZA, absRed_SZA ,label=["absorptance"])
plot(SZA, transRed_SZA,label=["transmittance"])
scatter(RAMI_SZA, RAMI_frefRed_050_SNW)
scatter(RAMI_SZA, RAMI_fabsRed_050_SNW)
scatter(RAMI_SZA, RAMI_ftranRed_050_SNW)
title("050_SNW - Default")
ylabel("Radiation Partitioning")
xlabel("Sun Zenith Angle")
xlim(0.0, 90.)
ylim(-0.05, 1.0)
xticks(0:20:91.)
yticks(0:0.2:1.0)

subplot(1,2,2)
plot(SZA, reflRed_clump_SZA,label=["reflectance"])
plot(SZA, absRed_clump_SZA ,label=["absorptance"])
plot(SZA,transRed_clump_SZA,label=["transmittance"])
scatter(RAMI_SZA, RAMI_frefRed_050_SNW,label="RAMI reflectance")
scatter(RAMI_SZA, RAMI_fabsRed_050_SNW,label="RAMI absorptance")
scatter(RAMI_SZA, RAMI_ftranRed_050_SNW,label="RAMI transmittance")
title("Clumping")
ylabel("Radiation Partitioning")
xlabel("Sun Zenith Angle")
xlim(0.0, 90.)
ylim(-0.05, 1.0)
xticks(0:20:91.)
yticks(0:0.2:1.0)

legend()
gcf()

##Defining all reference values for the Medium case

RAMI_SZA = [27.,60.,83.]


RAMI_fabsRed_150_BLK =  [0.28137804, 0.46514268999999997, 0.89063486]
RAMI_frefRed_150_BLK =  [0.00923676, 0.01379672, 0.02970703]
RAMI_ftranRed_150_BLK =  [0.7093851999999999, 0.52106059, 0.07965811]

RAMI_fabsRed_150_MED =  [0.31403827, 0.49003033, 0.89432051]
RAMI_frefRed_150_MED =  [0.06195053, 0.05151941, 0.03561715]
RAMI_ftranRed_150_MED =  [0.7104761399999999, 0.52197456, 0.07977039000000001]

RAMI_fabsRed_150_SNW =  [0.5431621799999999, 0.66519762, 0.9201217300000001]
RAMI_frefRed_150_SNW =  [0.43100610000000006, 0.31581022999999997, 0.07698033]
RAMI_ftranRed_150_SNW =  [0.71754777, 0.52755972, 0.08049832999999999]

##Medium case with black soil
reflRed_SZA,absRed_SZA,transRed_SZA,reflRed_clump_SZA,absRed_clump_SZA,transRed_clump_SZA=RAMI_case(1.5017, 0.0, 0.405417644);

SZA=0:1:85

figure(figsize=(10,5))

subplot(1,2,1)
plot(SZA, reflRed_SZA,label=["reflectance"])
plot(SZA, absRed_SZA ,label=["absorptance"])
plot(SZA, transRed_SZA,label=["transmittance"])
scatter(RAMI_SZA, RAMI_frefRed_150_BLK)
scatter(RAMI_SZA, RAMI_fabsRed_150_BLK)
scatter(RAMI_SZA, RAMI_ftranRed_150_BLK)
title("150_BLK - Default")
ylabel("Radiation Partitioning")
xlabel("Sun Zenith Angle")
xlim(0.0, 90.)
ylim(-0.05, 1.0)
xticks(0:20:91.)
yticks(0:0.2:1.0)

subplot(1,2,2)
plot(SZA, reflRed_clump_SZA,label=["reflectance"])
plot(SZA, absRed_clump_SZA ,label=["absorptance"])
plot(SZA,transRed_clump_SZA,label=["transmittance"])
scatter(RAMI_SZA, RAMI_frefRed_150_BLK,label="RAMI reflectance")
scatter(RAMI_SZA, RAMI_fabsRed_150_BLK,label="RAMI absorptance")
scatter(RAMI_SZA, RAMI_ftranRed_150_BLK,label="RAMI transmittance")
title("Clumping")
ylabel("Radiation Partitioning")
xlabel("Sun Zenith Angle")
xlim(0.0, 90.)
ylim(-0.05, 1.0)
xticks(0:20:91.)
yticks(0:0.2:1.0)

legend()
gcf()

##Medium case with medium soil
reflRed_SZA,absRed_SZA,transRed_SZA,reflRed_clump_SZA,absRed_clump_SZA,transRed_clump_SZA=RAMI_case(1.5017, 0.1217, 0.405417644);

SZA=0:1:85

figure(figsize=(10,5))

subplot(1,2,1)
plot(SZA, reflRed_SZA,label=["reflectance"])
plot(SZA, absRed_SZA ,label=["absorptance"])
plot(SZA, transRed_SZA,label=["transmittance"])
scatter(RAMI_SZA, RAMI_frefRed_150_MED)
scatter(RAMI_SZA, RAMI_fabsRed_150_MED)
scatter(RAMI_SZA, RAMI_ftranRed_150_MED)
title("150_MED - Default")
ylabel("Radiation Partitioning")
xlabel("Sun Zenith Angle")
xlim(0.0, 90.)
ylim(-0.05, 1.0)
xticks(0:20:91.)
yticks(0:0.2:1.0)

subplot(1,2,2)
plot(SZA, reflRed_clump_SZA,label=["reflectance"])
plot(SZA, absRed_clump_SZA ,label=["absorptance"])
plot(SZA,transRed_clump_SZA,label=["transmittance"])
scatter(RAMI_SZA, RAMI_frefRed_150_MED,label="RAMI reflectance")
scatter(RAMI_SZA, RAMI_fabsRed_150_MED,label="RAMI absorptance")
scatter(RAMI_SZA, RAMI_ftranRed_150_MED,label="RAMI transmittance")
title("Clumping")
ylabel("Radiation Partitioning")
xlabel("Sun Zenith Angle")
xlim(0.0, 90.)
ylim(-0.05, 1.0)
xticks(0:20:91.)
yticks(0:0.2:1.0)

legend()
gcf()

##Medium case with black soil
reflRed_SZA,absRed_SZA,transRed_SZA,reflRed_clump_SZA,absRed_clump_SZA,transRed_clump_SZA=RAMI_case(1.5017, 0.9640, 0.405417644);

SZA=0:1:85

figure(figsize=(10,5))

subplot(1,2,1)
plot(SZA, reflRed_SZA,label=["reflectance"])
plot(SZA, absRed_SZA ,label=["absorptance"])
plot(SZA, transRed_SZA,label=["transmittance"])
scatter(RAMI_SZA, RAMI_frefRed_150_SNW)
scatter(RAMI_SZA, RAMI_fabsRed_150_SNW)
scatter(RAMI_SZA, RAMI_ftranRed_150_SNW)
title("150_SNW - Default")
ylabel("Radiation Partitioning")
xlabel("Sun Zenith Angle")
xlim(0.0, 90.)
ylim(-0.05, 1.0)
xticks(0:20:91.)
yticks(0:0.2:1.0)

subplot(1,2,2)
plot(SZA, reflRed_clump_SZA,label=["reflectance"])
plot(SZA, absRed_clump_SZA ,label=["absorptance"])
plot(SZA,transRed_clump_SZA,label=["transmittance"])
scatter(RAMI_SZA, RAMI_frefRed_150_SNW,label="RAMI reflectance")
scatter(RAMI_SZA, RAMI_fabsRed_150_SNW,label="RAMI absorptance")
scatter(RAMI_SZA, RAMI_ftranRed_150_SNW,label="RAMI transmittance")
title("Clumping")
ylabel("Radiation Partitioning")
xlabel("Sun Zenith Angle")
xlim(0.0, 90.)
ylim(-0.05, 1.0)
xticks(0:20:91.)
yticks(0:0.2:1.0)

legend()
gcf()

##Defining all reference values for the Dense case

RAMI_SZA = [27.,60.,83.]


RAMI_fabsRed_250_BLK =  [0.46852539, 0.70426097, 0.9461774300000001]
RAMI_frefRed_250_BLK =  [0.01445858, 0.02016963, 0.03477486]
RAMI_ftranRed_250_BLK =  [0.51701603, 0.2755694, 0.01904771]

RAMI_fabsRed_250_MED =  [0.50540545, 0.72429659, 0.94742381]
RAMI_frefRed_250_MED =  [0.03953053, 0.03315039, 0.03580858]
RAMI_ftranRed_250_MED =  [0.51811911, 0.27616192, 0.01909098]

RAMI_fabsRed_250_SNW =  [0.76512258, 0.86538802, 0.9562473199999999]
RAMI_frefRed_250_SNW =  [0.21595537, 0.124503, 0.043056080000000004]
RAMI_ftranRed_250_SNW =  [0.5256125, 0.280805, 0.01935]

##Dense case with black soil
reflRed_SZA,absRed_SZA,transRed_SZA,reflRed_clump_SZA,absRed_clump_SZA,transRed_clump_SZA=RAMI_case(2.5007, 0.0, 0.45946608);

SZA=0:1:85

figure(figsize=(10,5))

subplot(1,2,1)
plot(SZA, reflRed_SZA,label=["reflectance"])
plot(SZA, absRed_SZA ,label=["absorptance"])
plot(SZA, transRed_SZA,label=["transmittance"])
scatter(RAMI_SZA, RAMI_frefRed_250_BLK)
scatter(RAMI_SZA, RAMI_fabsRed_250_BLK)
scatter(RAMI_SZA, RAMI_ftranRed_250_BLK)
title("250_BLK - Default")
ylabel("Radiation Partitioning")
xlabel("Sun Zenith Angle")
xlim(0.0, 90.)
ylim(-0.05, 1.0)
xticks(0:20:91.)
yticks(0:0.2:1.0)

subplot(1,2,2)
plot(SZA, reflRed_clump_SZA,label=["reflectance"])
plot(SZA, absRed_clump_SZA ,label=["absorptance"])
plot(SZA,transRed_clump_SZA,label=["transmittance"])
scatter(RAMI_SZA, RAMI_frefRed_250_BLK,label="RAMI reflectance")
scatter(RAMI_SZA, RAMI_fabsRed_250_BLK,label="RAMI absorptance")
scatter(RAMI_SZA, RAMI_ftranRed_250_BLK,label="RAMI transmittance")
title("Clumping")
ylabel("Radiation Partitioning")
xlabel("Sun Zenith Angle")
xlim(0.0, 90.)
ylim(-0.05, 1.0)
xticks(0:20:91.)
yticks(0:0.2:1.0)

legend()
gcf()

##Medium case with medium soil
reflRed_SZA,absRed_SZA,transRed_SZA,reflRed_clump_SZA,absRed_clump_SZA,transRed_clump_SZA=RAMI_case(2.5007, 0.1217, 0.45946608);

SZA=0:1:85

figure(figsize=(10,5))

subplot(1,2,1)
plot(SZA, reflRed_SZA,label=["reflectance"])
plot(SZA, absRed_SZA ,label=["absorptance"])
plot(SZA, transRed_SZA,label=["transmittance"])
scatter(RAMI_SZA, RAMI_frefRed_250_MED)
scatter(RAMI_SZA, RAMI_fabsRed_250_MED)
scatter(RAMI_SZA, RAMI_ftranRed_250_MED)
title("250_MED - Default")
ylabel("Radiation Partitioning")
xlabel("Sun Zenith Angle")
xlim(0.0, 90.)
ylim(-0.05, 1.0)
xticks(0:20:91.)
yticks(0:0.2:1.0)

subplot(1,2,2)
plot(SZA, reflRed_clump_SZA,label=["reflectance"])
plot(SZA, absRed_clump_SZA ,label=["absorptance"])
plot(SZA,transRed_clump_SZA,label=["transmittance"])
scatter(RAMI_SZA, RAMI_frefRed_250_MED,label="RAMI reflectance")
scatter(RAMI_SZA, RAMI_fabsRed_250_MED,label="RAMI absorptance")
scatter(RAMI_SZA, RAMI_ftranRed_250_MED,label="RAMI transmittance")
title("Clumping")
ylabel("Radiation Partitioning")
xlabel("Sun Zenith Angle")
xlim(0.0, 90.)
ylim(-0.05, 1.0)
xticks(0:20:91.)
yticks(0:0.2:1.0)

legend()
gcf()

##Dense case with snowy soil
reflRed_SZA,absRed_SZA,transRed_SZA,reflRed_clump_SZA,absRed_clump_SZA,transRed_clump_SZA=RAMI_case(2.5007, 0.9640, 0.45946608);

SZA=0:1:85

figure(figsize=(10,5))

subplot(1,2,1)
plot(SZA, reflRed_SZA,label=["reflectance"])
plot(SZA, absRed_SZA ,label=["absorptance"])
plot(SZA, transRed_SZA,label=["transmittance"])
scatter(RAMI_SZA, RAMI_frefRed_250_SNW)
scatter(RAMI_SZA, RAMI_fabsRed_250_SNW)
scatter(RAMI_SZA, RAMI_ftranRed_250_SNW)
title("250_SNW - Default")
ylabel("Radiation Partitioning")
xlabel("Sun Zenith Angle")
xlim(0.0, 90.)
ylim(-0.05, 1.0)
xticks(0:20:91.)
yticks(0:0.2:1.0)

subplot(1,2,2)
plot(SZA, reflRed_clump_SZA,label=["reflectance"])
plot(SZA, absRed_clump_SZA ,label=["absorptance"])
plot(SZA,transRed_clump_SZA,label=["transmittance"])
scatter(RAMI_SZA, RAMI_frefRed_250_SNW,label="RAMI reflectance")
scatter(RAMI_SZA, RAMI_fabsRed_250_SNW,label="RAMI absorptance")
scatter(RAMI_SZA, RAMI_ftranRed_250_SNW,label="RAMI transmittance")
title("Clumping")
ylabel("Radiation Partitioning")
xlabel("Sun Zenith Angle")
xlim(0.0, 90.)
ylim(-0.05, 1.0)
xticks(0:20:91.)
yticks(0:0.2:1.0)

legend()
gcf()

# Define a few wavelengths:
wl_blue = 450.0;
wl_red = 600.0;
wl_FarRed = 740.0;
wl_Red = 685.0;
ind_wle_blue  = argmin(abs.(wl_set.wle .-wl_blue));
ind_wle_red = argmin(abs.(wl_set.wle .-wl_red));
ind_wlf_FR  = argmin(abs.(wl_set.wlf .-wl_FarRed));
ind_wlf_R  = argmin(abs.(wl_set.wlf .-wl_Red));
ind_red = argmin(abs.(wl_set.wl .-wl_Red));
ind_NIR = argmin(abs.(wl_set.wl .-800));

SIF_FR = Float32[]
SIF_R = Float32[]
reflVIS = Float32[]
reflNIR = Float32[]

#### Just running the code over all geometries:
##MED
soil.albedo_SW[:] .=0.1217;
#### Set sun SZA to 27 degrees
angles.tts=27.
#### Set 0 azimuth (principal plane)
angles.psi=0

##Adding clumping
canopy_rt.Ω = 1.0
#### LAI of 3:
canopy_rt.LAI = 2.5007
#### Define VZA
VZA=collect(-89.5:0.5:89.5)

for VZA_ in VZA
    angles.tto=VZA_
    compute_canopy_geometry!(canopy_rt, angles, canOpt_rt)
    compute_canopy_matrices!(arrayOfLeaves, canOpt_rt);
    simulate_short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil);
    computeSIF_Fluxes!(arrayOfLeaves, canOpt_rt, canRad_rt, canopy_rt, soil, wl_set);
    #### Handpicked indices in
    push!(reflVIS, canRad_rt.alb_obs[ind_red])
    push!(reflNIR, canRad_rt.alb_obs[ind_NIR])
    push!(SIF_R , canRad_rt.SIF_obs[ind_wlf_R])
    push!(SIF_FR, canRad_rt.SIF_obs[ind_wlf_FR ])
end


##Adding clumping
canopy_rt.Ω = 0.45946608

SIF_FR_clump = Float32[]
SIF_R_clump = Float32[]
reflVIS_clump = Float32[]
reflNIR_clump = Float32[]


for VZA_ in VZA
    angles.tto=VZA_
    compute_canopy_geometry!(canopy_rt, angles, canOpt_rt)
    compute_canopy_matrices!(arrayOfLeaves, canOpt_rt);
    simulate_short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil);
    computeSIF_Fluxes!(arrayOfLeaves, canOpt_rt, canRad_rt, canopy_rt, soil, wl_set);
    #### Handpicked indices in
    push!(reflVIS_clump, canRad_rt.alb_obs[ind_red])
    push!(reflNIR_clump, canRad_rt.alb_obs[ind_NIR])
    push!(SIF_R_clump , canRad_rt.SIF_obs[ind_wlf_R])
    push!(SIF_FR_clump, canRad_rt.SIF_obs[ind_wlf_FR ])
end

#### Plots Visible
figure()
plot(VZA, reflVIS, color=:black,label="Red Reflectance", lw=2)
plot(VZA, SIF_R/30, color=:orange,label="Red SIF (/30)", lw=2)
plot(VZA, reflVIS_clump, color=:black, ls="--", lw=2,label="Red Reflectance w/ Clumping")
plot(VZA, SIF_R_clump/30, color=:orange, ls="--", lw=2,label="Red SIF (/30) w/ Clumping")
xlabel("Viewing Zenith Angle")
legend()
gcf()

#### Plots Visible
figure()
plot(VZA, reflNIR, color=:black,label="NIR Reflectance", lw=2)
plot(VZA, SIF_FR/6, color=:orange,label="Far Red SIF (/6)", lw=2)
plot(VZA, reflNIR_clump, color=:black, ls="--", lw=2,label="NIR Reflectance w/ Clumping")
plot(VZA, SIF_FR_clump/6, color=:orange, ls="--", lw=2,label="Far Red SIF (/6) w/ Clumping")
xlabel("Viewing Zenith Angle")
legend()
gcf()

reflVIS = Float32[]
reflNIR = Float32[]
SIF_FR = Float32[]
SIF_R  = Float32[]

##MED
soil.albedo_SW[:] .=0.1217;
angles.tts=27.
angles.psi=0
canopy_rt.LAI=2.5007
canopy_rt.Ω = 1.0
for psi=0:360
    angles.psi=psi
    for VZA=0:1:85
        angles.tto=VZA

        compute_canopy_geometry!(canopy_rt, angles, canOpt_rt)
        compute_canopy_matrices!(arrayOfLeaves, canOpt_rt);
        simulate_short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil);
        computeSIF_Fluxes!(arrayOfLeaves, canOpt_rt, canRad_rt, canopy_rt, soil, wl_set);
        push!(reflVIS, canRad_rt.alb_obs[28])
        push!(reflNIR, canRad_rt.alb_obs[52])
        push!(SIF_R , canRad_rt.SIF_obs[8])
        push!(SIF_FR, canRad_rt.SIF_obs[20])
    end
end

##Adding clumping
canopy_rt.Ω = 0.45946608

SIF_FR_clump = Float32[]
SIF_R_clump = Float32[]
reflVIS_clump = Float32[]
reflNIR_clump = Float32[]

for psi=0:360
    angles.psi=psi
    for VZA=0:1:85
        angles.tto=VZA

        compute_canopy_geometry!(canopy_rt, angles, canOpt_rt)
        compute_canopy_matrices!(arrayOfLeaves, canOpt_rt);
        simulate_short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil);
        computeSIF_Fluxes!(arrayOfLeaves, canOpt_rt, canRad_rt, canopy_rt, soil, wl_set);
        push!(reflVIS_clump, canRad_rt.alb_obs[28])
        push!(reflNIR_clump, canRad_rt.alb_obs[52])
        push!(SIF_R_clump , canRad_rt.SIF_obs[8])
        push!(SIF_FR_clump, canRad_rt.SIF_obs[20])
    end
end

A = reshape(reflNIR, ( 86,361));
B = reshape(reflVIS, ( 86,361));
SIFFER = reshape(SIF_R, ( 86,361));
SIFFER_FR = reshape(SIF_FR, ( 86,361));

A_clump = reshape(reflNIR_clump, ( 86,361));
B_clump = reshape(reflVIS_clump, ( 86,361));
SIFFER_clump = reshape(SIF_R_clump, ( 86,361));
SIFFER_FR_clump = reshape(SIF_FR_clump, ( 86,361));

figure(figsize=(10,5))
subplot(1,2,1, polar=true)
grid(false)
hm = contourf(deg2rad.(collect((0:360))),collect(0:1:85),  A,  cmap=:viridis, levels=collect(0.25:0.012:0.5))
title("NIR reflectance BRDF")
yticks([])
colorbar()
subplot(1,2,2, polar=true)
grid(false)
hm = contourf(deg2rad.(collect((0:360))),collect(0:1:85),  A_clump,  cmap=:viridis, levels=collect(0.25:0.012:0.5))
title("Clumping")
yticks([])
colorbar()
gcf()

figure(figsize=(10,5))
subplot(1,2,1, polar=true)
grid(false)
hm = contourf(deg2rad.(collect((0:360))),collect(0:1:85),  B,  cmap=:viridis, levels=collect(0.0:0.005:0.045))
title("Red reflectance BRDF")
yticks([])
colorbar()
subplot(1,2,2, polar=true)
grid(false)
hm = contourf(deg2rad.(collect((0:360))),collect(0:1:85),  B_clump,  cmap=:viridis, levels=collect(0.0:0.005:0.045))
title("Clumping")
yticks([])
colorbar()
gcf()

figure(figsize=(10,5))
subplot(1,2,1, polar=true)
grid(false)
hm = contourf(deg2rad.(collect((0:360))),collect(0:1:85),  SIFFER,  cmap=:viridis, levels=collect(0.3:0.05:0.8))
title("Red SIF emission BRDF")
yticks([])
colorbar()
subplot(1,2,2, polar=true)
grid(false)
hm = contourf(deg2rad.(collect((0:360))),collect(0:1:85),  SIFFER_clump,  cmap=:viridis, levels=collect(0.3:0.05:0.8))
title("Clumping")
yticks([])
colorbar()
gcf()

figure(figsize=(10,5))
subplot(1,2,1, polar=true)
grid(false)
hm = contourf(deg2rad.(collect((0:360))),collect(0:1:85),  SIFFER_FR,  cmap=:viridis, levels=collect(1.:0.1:3.5))
title("Far-Red SIF emission BRDF")
yticks([])
colorbar()
subplot(1,2,2, polar=true)
grid(false)
hm = contourf(deg2rad.(collect((0:360))),collect(0:1:85),  SIFFER_FR_clump,  cmap=:viridis, levels=collect(1.:0.1:3.5))
title("Clumping")
yticks([])
colorbar()
gcf()

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

