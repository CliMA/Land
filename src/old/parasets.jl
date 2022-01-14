###############################################################################
#
# Temperature Dedenpdency Parameter sets
# Data source: Bernacchi et al. (2001)
# Improved temperature response functions for models of Rubisco‐limited
# photosynthesis
#
###############################################################################
""" [`ArrheniusTD`](@ref) type Kc TD from Bernacchi's data """
KcTDBernacchi(FT) = Arrhenius{FT}(T_25(), 41.0264925, 79430.0);

""" [`ArrheniusTD`](@ref) type Ko TD from Bernacchi's data """
KoTDBernacchi(FT) = Arrhenius{FT}(T_25(), 28208.88, 36380.0);

""" [`ArrheniusTD`](@ref) type Respiration TD from Bernacchi's data """
RespirationTDBernacchi(FT) = Arrhenius{FT}(T_25(), 1, 46390.0);

""" [`ArrheniusTD`](@ref) type Vcmax TD from Bernacchi's data """
VcmaxTDBernacchi(FT) = Arrhenius{FT}(T_25(), 1, 65330.0);

""" [`ArrheniusTD`](@ref) type Vomax TD from Bernacchi's data """
VomaxTDBernacchi(FT) = Arrhenius{FT}(T_25(), 1, 60110.0);

""" [`ArrheniusTD`](@ref) type ``Γ^{*}`` TD from Bernacchi's data """
ΓStarTDBernacchi(FT) = Arrhenius{FT}(T_25(), 4.33164375, 37830.0);


###############################################################################
#
# Temperature Dedenpdency Parameter sets
# Data source: Boyd et al. (2001)
# Temperature responses of C4 photosynthesis: biochemical analysis of Rubisco,
# phosphoenolpyruvate carboxylase, and carbonic anhydrase in Setaria viridis
#
###############################################################################
""" [`ArrheniusTD`](@ref) type Kpep TD from Boyd's data """
KpepTDBoyd(FT) = Arrhenius{FT}(T_25(), 16.0, 36300.0);

""" [`ArrheniusPeakTD`](@ref) type Vpmax TD from Boyd's data """
VpmaxTDBoyd(FT) = ArrheniusPeak{FT}(T_25(), 1, 94800.0, 73300.0, 250.0);


###############################################################################
#
# Temperature Dedenpdency Parameter sets
# Data source: Leuning (2002)
# Temperature dependence of two parameters in a photosynthesis model
#
###############################################################################
""" [`ArrheniusPeakTD`](@ref) type Jmax TD from Leuning's data """
JmaxTDLeuning(FT) =  ArrheniusPeak{FT}(T_25(), 1, 50300.0, 152044.0, 495.0);

""" [`ArrheniusPeakTD`](@ref) type Vcmax TD from Leuning's data """
VcmaxTDLeuning(FT) = ArrheniusPeak{FT}(T_25(), 1, 73637.0, 149252.0, 486.0);


###############################################################################
#
# Temperature Dedenpdency Parameter sets
# Data source: Lavigne and Ryan (1997)
# Growth and maintenance respiration rates of aspen, blackspruce and jack pine
#     stems at northern and southern BOREAS sites
#
###############################################################################
""" [`Q10TD`](@ref) type Respiration TD for angiosperms per biomass """
Q10TDAngiosperm(FT) = Q10{FT}(T_25(), 0.014/8760, 1.4)

""" [`Q10TD`](@ref) type Respiration TD for symnosperms per biomass """
Q10TDGymnosperm(FT) = Q10{FT}(T_25(), 0.0425/8760, 1.7)


###############################################################################
#
# Temperature Dedenpdency Parameter sets
# Data source missing
#
###############################################################################
""" [`ArrheniusTD`](@ref) type Kc TD """
KcTDCLM(FT) = Arrhenius{FT}(T_25(), 40.49, 79430.0);

""" [`ArrheniusTD`](@ref) type Ko TD """
KoTDCLM(FT) = Arrhenius{FT}(T_25(), 27840.0, 36380.0);

""" [`ArrheniusTD`](@ref) type Kpep TD """
KpepTDCLM(FT) = Arrhenius{FT}(T_25(), 8.0, 36000.0);

""" [`ArrheniusTD`](@ref) type Γ* TD """
ΓStarTDCLM(FT) = Arrhenius{FT}(T_25(), 4.275, 37830.0);

""" [`ArrheniusPeakTD`](@ref) type Jmax TD """
JmaxTDBernacchi(FT) = ArrheniusPeak{FT}(T_25(), 1, 57500.0, 439000.0, 1400.0);

""" [`ArrheniusPeakTD`](@ref) type Jmax TD """
JmaxTDCLM(FT) = ArrheniusPeak{FT}(T_25(), 1, 43540.0, 150000.0, 490.0);

""" [`ArrheniusPeakTD`](@ref) type Respiration TD """
RespirationTDCLM(FT) = ArrheniusPeak{FT}(T_25(), 1, 46390.0, 150650.0, 490.0);

""" [`ArrheniusPeakTD`](@ref) type Vcmax TD """
VcmaxTDCLM(FT) = ArrheniusPeak{FT}(T_25(), 1, 65330.0, 150000.0, 490.0);


###############################################################################
#
# Vcmax to R correlation
#
###############################################################################
""" A constant of 0.01 """
VtoRCollatz(FT) = FT(0.010);

""" A constant of 0.015 """
VtoRDefault(FT) = FT(0.015);


###############################################################################
#
# Fluorescence model parameter set
# Data source: van der Tol et al. (2014)
# Models of fluorescence and photosynthesis for interpreting measurements of
#     solar-induced chlorophyll fluorescence
#
###############################################################################
""" [`FluoParaSet`](@ref) type parameter set using all data """
FluorescenceVanDerTol(FT) = FluoParaSet{FT}(2.48, 2.83, 0.114);

""" [`FluoParaSet`](@ref) type parameter set using Flexas's data (drought) """
FluorescenceVanDerTolDrought(FT) = FluoParaSet{FT}(5.01, 1.93, 10.0);


###############################################################################
#
# Pre-set photosynthesis model parameter sets
#
###############################################################################
""" [`C3ParaSet`](@ref) type C3 photosynthesis using Bernacchi's data """
function C3Bernacchi(FT)
    JT  = JmaxTDBernacchi(FT);
    KcT = KcTDBernacchi(FT);
    KoT = KoTDBernacchi(FT);
    ReT = RespirationTDBernacchi(FT);
    VcT = VcmaxTDBernacchi(FT);
    ΓsT = ΓStarTDBernacchi(FT);
    Flu = FluorescenceVanDerTolDrought(FT);
    VR  = VtoRDefault(FT);
    E1  = FT(4);
    E2  = FT(8);
    return C3ParaSet{FT}(JT, KcT, KoT, ReT, VcT, ΓsT, Flu, VR, E1, E2)
end

""" [`C3ParaSet`](@ref) type C3 photosynthesis using CLM5's data """
function C3CLM(FT)
    JT  = JmaxTDCLM(FT);
    KcT = KcTDCLM(FT);
    KoT = KoTDCLM(FT);
    ReT = RespirationTDCLM(FT);
    VcT = VcmaxTDCLM(FT);
    ΓsT = ΓStarTDCLM(FT);
    Flu = FluorescenceVanDerTolDrought(FT);
    VR  = VtoRDefault(FT);
    E1  = FT(4);
    E2  = FT(8);
    return C3ParaSet{FT}(JT, KcT, KoT, ReT, VcT, ΓsT, Flu, VR, E1, E2)
end

""" [`C4ParaSet`](@ref) type C4 photosynthesis using CLM5's data """
function C4CLM(FT)
    KpT = KpepTDCLM(FT);
    ReT = RespirationTDCLM(FT);
    VcT = VcmaxTDCLM(FT);
    VpT = VpmaxTDBoyd(FT);
    Flu = FluorescenceVanDerTolDrought(FT);
    VR  = VtoRDefault(FT);
    return C4ParaSet{FT}(KpT, ReT, VcT, VpT, Flu, VR)
end
