# Load land.Photosynthesis module
using Land.Hydraulics
using Land.Photosynthesis

# Load other packages
using PyPlot

# Use Float32 by default, if Float32 runs OK, then no problem with Float64
FT = Float32;




# Background of photosynthesis
# Photosynthesis types: C3/C4/CAM
#
# Differences of the photosynthesis types
# C3:  Light capture and CO2 fixation in the same cell
#            happen simutaneously
# C4:  Light capture in bundle sheath (with chloroplasts)
#            CO2 fixation in mesophyll (no chloroplast)
#            happen simultaneously
# CAM: Light capture and CO2 fixation in the same cell
#            Light capture in the day
#            CO2 fixation at night
# Some plants are able to do
#            C3  metabolism when water is sufficient
#            CAM metabolism when water is inadequate
#            e.g., orchids
#
# Tips of C3/C4 photosynthesis
# * When C3 photosynthesis originated, CO2 concentration was super high,
#            but O2 concentration was very low
#            "plants" did not the RubisCO can fix O2 as well
#            This is what the "C" and "O" mean
#            Ribulose-1,5-bisphosphate carboxylase/oxygenase
# * C3 plants has PEPcase as well
#            but not working like in C4 plants
# * With the decreasing of [CO2] but increasing of [O2]
#            O2 fixation compete with CO2 fixation
#            This is called photorespiration (about 1/3 CO2 fixation)
#            A threshold where O2 and CO2 fixation rates equal
# * There is a driving force to photosynthesize at low [CO2]
#            Solution 1: Evolve new enzyme to fix CO2 only
#            Solution 2: Reduce [O2] in chloroplast
#            Solution 3: Increase [CO2 ] in chloroplasts
#            S1 is too expensive
#            S2 is not phyiologically possible because CO2 > O2 > H2O
#            S3 PEPcase already there, only minor structure changes requried
# * At the end, plants choose to do a minor revision
#            Revision 1: separate metabolism in space --> C4  metabolism
#            Revision 2: separate metabolism in time  --> CAM metabolism




# C3 photosynthesis
# Dark reaction limited --- CO2 and O2 competition for RubisCO
# 
#                      Pi - Γ*
# Ac = Vcmax * -----------------------
#                Pi + Kc * (1+O2/Ko)
#
# Light reaction limited --- electron transport
#
#             Pi - Γ*
# Aj = J * -------------
#            4Pi + 8Γ*
#
# J = f(Jmax, PAR, Quantum Yield)
#
# Product limited
#
# Ap = a * Vcmax
#
# Final rates
# 
# Ag = min(Ac, Aj, Ap)
# 
# An = Ag - Rd
#



#=
# All reactions are impacted by temperature
# https://clima.github.io/Land/dev/pages/Photosynthesis/#Temperature-Dependency-1

_TD      = Photosynthesis.RespirationTDBernacchi(FT);
_TD_peak = Photosynthesis.RespirationTDCLM(FT);
@show typeof(_TD);
@show typeof(_TD_peak);

_Tl = collect(FT, 1:40) .+ FT(273.15);
_C  = Photosynthesis.arrhenius_correction.([_TD]     , _Tl);
_Cp = Photosynthesis.arrhenius_correction.([_TD_peak], _Tl);

figure();
tight_layout(true);
plot(_Tl .- FT(273.15), _C , "r-", label="TD"     );
plot(_Tl .- FT(273.15), _Cp, "b-", label="Peak TD");
xlabel("\$T_{leaf}\$ (\$^{\\circ}\$C)"      , fontsize=16);
ylabel("Value relative to 25 \$^{\\circ}\$C", fontsize=16);
legend();
gcf();




# CO2 impact on photosynthetic rates
# Use public functions from Photosynthesis module
# create standard leaves for C3/C4 photosynthesis
# https://clima.github.io/Land/dev/pages/Photosynthesis/#Leaf-and-Environment-1
leaf_3 = Leaf{FT}(APAR=1000);
leaf_4 = Leaf{FT}(APAR=1000);

# create standard environmental conditions
envir = AirLayer{FT}();

# define the photosynthesis parasets
# https://clima.github.io/Land/dev/pages/Photosynthesis/#Parameter-Sets-Collection-1
c3_set = C3CLM(FT);
c4_set = C4CLM(FT);
c4_set.KpT = Photosynthesis.ArrheniusTD{FT}(16.0, 4329.806249108331, 14.522241318491803)

# Run photosynthesis model for different [CO2]
_Cl = collect(FT,1:200);

leaf_3.T = 310;
leaf_4.T = 310;

Ac_C3 = similar(_Cl); Aj_C3 = similar(_Cl);
Ap_C3 = similar(_Cl); Ag_C3 = similar(_Cl);
Ac_C4 = similar(_Cl); Aj_C4 = similar(_Cl);
Ap_C4 = similar(_Cl); Ag_C4 = similar(_Cl);

for i in eachindex(_Cl)
    leaf_3.p_i = _Cl[i];
    leaf_4.p_i = _Cl[i];
    leaf_photo_from_pi!(c3_set, leaf_3, envir);
    leaf_photo_from_pi!(c4_set, leaf_4, envir);
    
    Ac_C3[i] = leaf_3.Ac; Aj_C3[i] = leaf_3.Aj;
    Ap_C3[i] = leaf_3.Ap; Ag_C3[i] = leaf_3.Ag;
    Ac_C4[i] = leaf_4.Ac; Aj_C4[i] = leaf_4.Aj;
    Ap_C4[i] = leaf_4.Ap; Ag_C4[i] = leaf_4.Ag;
end

# plot the results
figure(figsize=(10,5));
tight_layout(true);

subplot(1,2,1);
plot(_Cl, Ac_C3, "g-", lw=2, label="Ac");
plot(_Cl, Aj_C3, "r-", lw=2, label="Aj");
plot(_Cl, Ap_C3, "b-", lw=2, label="Ap");
plot(_Cl, Ag_C3, "k-", lw=2, label="Ag");
xlabel("Leaf Internal CO₂ (Pa)", fontsize=16);
ylabel("A (μmol m⁻² s⁻¹)"      , fontsize=16);
legend(loc="lower right");

subplot(1,2,2);
plot(_Cl, Ac_C4, "g-", lw=2, label="Ac");
plot(_Cl, Aj_C4, "r-", lw=2, label="Aj");
plot(_Cl, Ap_C4, "b-", lw=2, label="Ap");
plot(_Cl, Ag_C4, "k-", lw=2, label="Ag");
xlabel("Leaf Internal CO₂ (Pa)", fontsize=16);
legend(loc="lower right");

gcf();




# PAR impact on photosynthetic rates
# Use public functions from Photosynthesis module

# create standard leaves for C3/C4 photosynthesis
leaf_3 = Leaf{FT}(p_i=20);
leaf_4 = Leaf{FT}(p_i=20);

# create standard environmental conditions
envir = AirLayer{FT}();

# define the photosynthesis parasets
c3_set = C3CLM(FT);
c4_set = C4CLM(FT);

# initialize the TD and Radiation
leaf_temperature_dependence!(c3_set, leaf_3, envir);
leaf_temperature_dependence!(c4_set, leaf_4, envir);

# Run photosynthesis model for different [CO2]
_Pl = collect(FT,0:10:2000);

Ac_C3 = similar(_Pl); Aj_C3 = similar(_Pl);
Ap_C3 = similar(_Pl); Ag_C3 = similar(_Pl);
Ac_C4 = similar(_Pl); Aj_C4 = similar(_Pl);
Ap_C4 = similar(_Pl); Ag_C4 = similar(_Pl);

for i in eachindex(_Pl)
    leaf_3.APAR = _Pl[i];
    leaf_4.APAR = _Pl[i];
    leaf_photo_from_pi!(c3_set, leaf_3, envir);
    leaf_photo_from_pi!(c4_set, leaf_4, envir);
    
    Ac_C3[i] = leaf_3.Ac; Aj_C3[i] = leaf_3.Aj;
    Ap_C3[i] = leaf_3.Ap; Ag_C3[i] = leaf_3.Ag;
    Ac_C4[i] = leaf_4.Ac; Aj_C4[i] = leaf_4.Aj;
    Ap_C4[i] = leaf_4.Ap; Ag_C4[i] = leaf_4.Ag;
end

# plot the results
figure(figsize=(10,5));
tight_layout(true);

subplot(1,2,1);
plot(_Pl, Ac_C3, "g-", lw=2, label="Ac");
plot(_Pl, Aj_C3, "r-", lw=2, label="Aj");
plot(_Pl, Ap_C3, "b-", lw=2, label="Ap");
plot(_Pl, Ag_C3, "k-", lw=2, label="Ag");
xlabel("Absorbed PAR (μmol m⁻² s⁻¹)", fontsize=16);
ylabel("A (μmol m⁻² s⁻¹)"           , fontsize=16);
legend(loc="lower right");

subplot(1,2,2);
plot(_Pl, Ac_C4, "g-", lw=2, label="Ac");
plot(_Pl, Aj_C4, "r-", lw=2, label="Aj");
plot(_Pl, Ap_C4, "b-", lw=2, label="Ap");
plot(_Pl, Ag_C4, "k-", lw=2, label="Ag");
xlabel("Absorbed PAR (μmol m⁻² s⁻¹)", fontsize=16);
legend(loc="lower right");

gcf()




# Stomatal responses to the environment, e.g., PAR
# Use public functions from Photosynthesis module
# https://clima.github.io/Land/dev/pages/Photosynthesis/#Stomatal-Response-to-Environment-1

# create standard leaves for C3/C4 photosynthesis
leaf_3 = Leaf{FT}();
leaf_4 = Leaf{FT}();

# create standard environmental conditions
envir = AirLayer{FT}();

# define the photosynthesis parasets
c3_set = C3CLM(FT);
c4_set = C4CLM(FT);
c3_set.Sto = Photosynthesis.ESMBallBerry{FT}(g1 = 16);
c4_set.Sto = Photosynthesis.ESMBallBerry{FT}(g1 = 8);

# Run photosynthesis model for different APAR
_Pl = collect(FT,0:10:1000);

Ac_C3 = similar(_Pl); Aj_C3 = similar(_Pl);
Ap_C3 = similar(_Pl); Ag_C3 = similar(_Pl);
An_C3 = similar(_Pl); gs_C3 = similar(_Pl);
Ac_C4 = similar(_Pl); Aj_C4 = similar(_Pl);
Ap_C4 = similar(_Pl); Ag_C4 = similar(_Pl);
An_C4 = similar(_Pl); gs_C4 = similar(_Pl);

for i in eachindex(_Pl)
    leaf_3.APAR = _Pl[i];
    leaf_4.APAR = _Pl[i];
    leaf_photo_from_envir!(c3_set, leaf_3, envir, c3_set.Sto);
    leaf_photo_from_envir!(c4_set, leaf_4, envir, c4_set.Sto);
    
    Ac_C3[i] = leaf_3.Ac; Aj_C3[i] = leaf_3.Aj  ;
    Ap_C3[i] = leaf_3.Ap; Ag_C3[i] = leaf_3.Ag  ;
    An_C3[i] = leaf_3.An; gs_C3[i] = leaf_3.g_sw;
    Ac_C4[i] = leaf_4.Ac; Aj_C4[i] = leaf_4.Aj  ;
    Ap_C4[i] = leaf_4.Ap; Ag_C4[i] = leaf_4.Ag  ;
    An_C4[i] = leaf_4.An; gs_C4[i] = leaf_4.g_sw;
end

# plot the results
figure(figsize=(15,10));
tight_layout(true);

subplot(2,2,1);
plot(_Pl, Ac_C3, "g-", lw=2, label="Ac");
plot(_Pl, Aj_C3, "r-", lw=2, label="Aj");
plot(_Pl, Ap_C3, "b-", lw=2, label="Ap");
plot(_Pl, Ag_C3, "k-", lw=2, label="Ag");
ylabel("A (μmol m⁻² s⁻¹)"           , fontsize=16);
legend(loc="lower right");

subplot(2,2,3);
plot(_Pl, Ac_C4, "g-", lw=2, label="Ac");
plot(_Pl, Aj_C4, "r-", lw=2, label="Aj");
plot(_Pl, Ap_C4, "b-", lw=2, label="Ap");
plot(_Pl, Ag_C4, "k-", lw=2, label="Ag");
xlabel("Absorbed PAR (μmol m⁻² s⁻¹)", fontsize=16);
ylabel("A (μmol m⁻² s⁻¹)"           , fontsize=16);
legend(loc="lower right");

subplot(1,2,2);
plot(_Pl, gs_C3, "k-", lw=2; label="C3");
plot(_Pl, gs_C4, "k:", lw=2;label="C4");
xlabel("Absorbed PAR (μmol m⁻² s⁻¹)", fontsize=16);
ylabel("Gsw (mol m⁻² s⁻¹)"          , fontsize=16);
legend(loc="lower right");

gcf();
=#



# Use another stomatal model scheme, e.g., optimization
# Empirical models: ESMBallBerry, ESMGentine, ESMLeuning, ESMMedlyn
# Optimization models: OSMEller, OSMSperry, OSMWang, OSMWAP

# create standard leaves for C3/C4 photosynthesis
leaf_3 = Leaf{FT}(hs=LeafHydraulics{FT}(k_sla=0.05));
leaf_4 = Leaf{FT}(hs=LeafHydraulics{FT}(k_sla=0.05));

# create standard environmental conditions
envir = AirLayer{FT}();

c3_set = C3CLM(FT);
c4_set = C4CLM(FT);
c3_set.Sto = Photosynthesis.OSMEller();
c4_set.Sto = Photosynthesis.OSMEller();

# Run photosynthesis model for different APAR
_Pl = collect(FT,0:10:1000);

Ac_C3 = similar(_Pl); Aj_C3 = similar(_Pl);
Ap_C3 = similar(_Pl); Ag_C3 = similar(_Pl);
An_C3 = similar(_Pl); gs_C3 = similar(_Pl);
Ac_C4 = similar(_Pl); Aj_C4 = similar(_Pl);
Ap_C4 = similar(_Pl); Ag_C4 = similar(_Pl);
An_C4 = similar(_Pl); gs_C4 = similar(_Pl);

# Run photosynthesis model for different APAR

for i in eachindex(_Pl)
    leaf_3.APAR = _Pl[i];
    leaf_4.APAR = _Pl[i];
    leaf_photo_from_envir!(c3_set, leaf_3, envir, c3_set.Sto);
    leaf_photo_from_envir!(c4_set, leaf_4, envir, c4_set.Sto);
    
    Ac_C3[i] = leaf_3.Ac; Aj_C3[i] = leaf_3.Aj  ;
    Ap_C3[i] = leaf_3.Ap; Ag_C3[i] = leaf_3.Ag  ;
    An_C3[i] = leaf_3.An; gs_C3[i] = leaf_3.g_sw;
    Ac_C4[i] = leaf_4.Ac; Aj_C4[i] = leaf_4.Aj  ;
    Ap_C4[i] = leaf_4.Ap; Ag_C4[i] = leaf_4.Ag  ;
    An_C4[i] = leaf_4.An; gs_C4[i] = leaf_4.g_sw;
end

# plot the results
figure(figsize=(15,10));
tight_layout(true);

subplot(2,2,1);
plot(_Pl, Ac_C3, "g-", lw=2, label="Ac");
plot(_Pl, Aj_C3, "r-", lw=2, label="Aj");
plot(_Pl, Ap_C3, "b-", lw=2, label="Ap");
plot(_Pl, Ag_C3, "k-", lw=2, label="Ag");
ylabel("A (μmol m⁻² s⁻¹)"           , fontsize=16);
legend(loc="lower right");

subplot(2,2,3);
plot(_Pl, Ac_C4, "g-", lw=2, label="Ac");
plot(_Pl, Aj_C4, "r-", lw=2, label="Aj");
plot(_Pl, Ap_C4, "b-", lw=2, label="Ap");
plot(_Pl, Ag_C4, "k-", lw=2, label="Ag");
xlabel("Absorbed PAR (μmol m⁻² s⁻¹)", fontsize=16);
ylabel("A (μmol m⁻² s⁻¹)"           , fontsize=16);
legend(loc="lower right");

subplot(1,2,2);
plot(_Pl, gs_C3, "k-", lw=2; label="C3");
plot(_Pl, gs_C4, "k:", lw=2;label="C4");
xlabel("Absorbed PAR (μmol m⁻² s⁻¹)", fontsize=16);
ylabel("Gsw (mol m⁻² s⁻¹)"          , fontsize=16);
legend(loc="lower right");

gcf();




# Response to drought
# create standard leaves for C3/C4 photosynthesis
#c3_set.Sto = Photosynthesis.ESMBallBerry{FT}(g1 = 16);
#c4_set.Sto = Photosynthesis.ESMBallBerry{FT}(g1 = 8);
#c3_set.Sto = Photosynthesis.ESMLeuning{FT}(g1 = 16);
#c4_set.Sto = Photosynthesis.ESMLeuning{FT}(g1 = 8);
#c3_set.Sto = Photosynthesis.ESMMedlyn{FT}(g1 = 16);
#c4_set.Sto = Photosynthesis.ESMMedlyn{FT}(g1 = 8);
#c3_set.Sto = Photosynthesis.ESMGentine{FT}(g1=16);
#c4_set.Sto = Photosynthesis.ESMGentine{FT}(g1=8);
#c3_set.Sto = Photosynthesis.OSMEller();
#c4_set.Sto = Photosynthesis.OSMEller();
#c3_set.Sto = Photosynthesis.OSMSperry();
#c4_set.Sto = Photosynthesis.OSMSperry();
#c3_set.Sto = Photosynthesis.OSMWang();
#c4_set.Sto = Photosynthesis.OSMWang();
# OSMWAP model and OMSWAPMod not tested yet
#c3_set.Sto = Photosynthesis.OSMWAP{FT}();
#c4_set.Sto = Photosynthesis.OSMWAP{FT}();

# simulate soil drought
_Sl = collect(FT,0:-0.01:leaf_3.hs.p_crt);

Ac_C3 = similar(_Sl); Aj_C3 = similar(_Sl);
Ap_C3 = similar(_Sl); Ag_C3 = similar(_Sl);
An_C3 = similar(_Sl); gs_C3 = similar(_Sl);
ec_C3 = similar(_Sl);
Ac_C4 = similar(_Sl); Aj_C4 = similar(_Sl);
Ap_C4 = similar(_Sl); Ag_C4 = similar(_Sl);
An_C4 = similar(_Sl); gs_C4 = similar(_Sl);
ec_C4 = similar(_Sl);

for i in eachindex(_Sl)
    leaf_3.p_ups = _Sl[i];
    leaf_4.p_ups = _Sl[i];
    leaf_3.hs.p_ups = _Sl[i];
    leaf_4.hs.p_ups = _Sl[i];
    leaf_photo_from_envir!(c3_set, leaf_3, envir, c3_set.Sto);
    leaf_photo_from_envir!(c4_set, leaf_4, envir, c4_set.Sto);
    
    Ac_C3[i] = leaf_3.Ac; Aj_C3[i] = leaf_3.Aj  ;
    Ap_C3[i] = leaf_3.Ap; Ag_C3[i] = leaf_3.Ag  ;
    An_C3[i] = leaf_3.An; gs_C3[i] = leaf_3.g_sw;
    ec_C3[i] = leaf_3.ec;
    Ac_C4[i] = leaf_4.Ac; Aj_C4[i] = leaf_4.Aj  ;
    Ap_C4[i] = leaf_4.Ap; Ag_C4[i] = leaf_4.Ag  ;
    An_C4[i] = leaf_4.An; gs_C4[i] = leaf_4.g_sw;
    ec_C4[i] = leaf_4.ec;
end

# plot the results
figure(figsize=(16,8));
tight_layout(true);

subplot(2,2,1);
plot(_Sl, Ac_C3, "g-", lw=2, label="Ac");
plot(_Sl, Aj_C3, "r-", lw=2, label="Aj");
plot(_Sl, Ap_C3, "b-", lw=2, label="Ap");
plot(_Sl, Ag_C3, "k-", lw=2, label="Ag");
ylabel("A (μmol m⁻² s⁻¹)"          , fontsize=16);
legend(loc="lower right");

subplot(2,2,3);
plot(_Sl, Ac_C4, "g-", lw=2, label="Ac");
plot(_Sl, Aj_C4, "r-", lw=2, label="Aj");
plot(_Sl, Ap_C4, "b-", lw=2, label="Ap");
plot(_Sl, Ag_C4, "k-", lw=2, label="Ag");
xlabel("Soil water potential (MPa)", fontsize=16);
ylabel("A (μmol m⁻² s⁻¹)"          , fontsize=16);
legend(loc="lower right");

subplot(1,2,2);
plot(_Sl, gs_C3, "k-", lw=2; label="C3");
plot(_Sl, gs_C4, "k:", lw=2;label="C4");
xlabel("Soil water potential (MPa)", fontsize=16);
ylabel("Gsw (mol m⁻² s⁻¹)"         , fontsize=16);
legend(loc="lower right");

gcf();

