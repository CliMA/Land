###############################################################################
#
# Temperature Dedenpdency Parameter sets
# Data source: Bernacchi et al. (2001)
# Improved temperature response functions for models of Rubisco‐limited
# photosynthesis
#
###############################################################################
""" [`ArrheniusTD`](@ref) type Kc TD from Bernacchi's data """
KcTDBernacchi(FT) =
    ArrheniusTD{FT}(41.0264925, 79430.0 / GAS_R(), 79430.0 / RT_25(FT));

""" [`ArrheniusTD`](@ref) type Ko TD from Bernacchi's data """
KoTDBernacchi(FT) =
    ArrheniusTD{FT}(  28208.88, 36380.0 / GAS_R(), 36380.0 / RT_25(FT));

""" [`ArrheniusTD`](@ref) type Respiration TD from Bernacchi's data """
RespirationTDBernacchi(FT) =
    ArrheniusTD{FT}(       Inf, 46390.0 / GAS_R(), 46390.0 / RT_25(FT));

""" [`ArrheniusTD`](@ref) type Vcmax TD from Bernacchi's data """
VcmaxTDBernacchi(FT) =
    ArrheniusTD{FT}(       Inf, 65330.0 / GAS_R(), 65330.0 / RT_25(FT));

""" [`ArrheniusTD`](@ref) type Vomax TD from Bernacchi's data """
VomaxTDBernacchi(FT) =
    ArrheniusTD{FT}(       Inf, 60110.0 / GAS_R(), 60110.0 / RT_25(FT));

""" [`ArrheniusTD`](@ref) type ``Γ^{*}`` TD from Bernacchi's data """
ΓStarTDBernacchi(FT) =
    ArrheniusTD{FT}(4.33164375, 37830.0 / GAS_R(), 37830.0 / RT_25(FT));








###############################################################################
#
# Temperature Dedenpdency Parameter sets
# Data source: Boyd et al. (2001)
# Temperature responses of C4 photosynthesis: biochemical analysis of Rubisco,
# phosphoenolpyruvate carboxylase, and carbonic anhydrase in Setaria viridis
#
###############################################################################
""" [`ArrheniusTD`](@ref) type Kpep TD from Boyd's data """
KpepTDBoyd(FT) =
    ArrheniusTD{FT}(16.0, 36300.0 / GAS_R(), 36300.0 / RT_25(FT));

""" [`ArrheniusPeakTD`](@ref) type Vpmax TD from Boyd's data """
function VpmaxTDBoyd(FT)
    ΔHa_to_RT25::FT = 94800.0 / RT_25(FT);
    ΔHd_to_R::FT = 73300.0 / GAS_R();
    ΔSv_to_R::FT = 250.0   / GAS_R();
    C::FT = 1 + exp( ΔSv_to_R - ΔHd_to_R/T_25(FT) );
    return ArrheniusPeakTD{FT}(ΔHa_to_RT25, ΔHd_to_R, ΔSv_to_R, C)
end








###############################################################################
#
# Temperature Dedenpdency Parameter sets
# Data source: Leuning (2002)
# Temperature dependence of two parameters in a photosynthesis model
#
###############################################################################
""" [`ArrheniusPeakTD`](@ref) type Jmax TD from Leuning's data """
function JmaxTDLeuning(FT)
    ΔHa_to_RT25::FT = 50300.0  / RT_25(FT);
    ΔHd_to_R::FT = 152044.0 / GAS_R();
    ΔSv_to_R::FT = 495.0    / GAS_R();
    C::FT = 1 + exp( ΔSv_to_R - ΔHd_to_R/T_25(FT) );
    return ArrheniusPeakTD{FT}(ΔHa_to_RT25, ΔHd_to_R, ΔSv_to_R, C)
end

""" [`ArrheniusPeakTD`](@ref) type Vcmax TD from Leuning's data """
function VcmaxTDLeuning(FT)
    ΔHa_to_RT25::FT = 73637.0  / RT_25(FT);
    ΔHd_to_R::FT = 149252.0 / GAS_R();
    ΔSv_to_R::FT = 486.0    / GAS_R();
    C::FT = 1 + exp( ΔSv_to_R - ΔHd_to_R/T_25(FT) );
    return ArrheniusPeakTD{FT}(ΔHa_to_RT25, ΔHd_to_R, ΔSv_to_R, C)
end








###############################################################################
#
# Temperature Dedenpdency Parameter sets
# Data source: Lavigne and Ryan (1997)
# Growth and maintenance respiration rates of aspen, blackspruce and jack pine
#     stems at northern and southern BOREAS sites
#
###############################################################################
""" [`Q10TD`](@ref) type Respiration TD for angiosperms per biomass """
Q10TDAngiosperm(FT) = Q10TD{FT}(0.014/8760, 298.15, 1.4)

""" [`Q10TD`](@ref) type Respiration TD for symnosperms per biomass """
Q10TDGymnosperm(FT) = Q10TD{FT}(0.0425/8760, 298.15, 1.7)








###############################################################################
#
# Temperature Dedenpdency Parameter sets
# Data source missing
#
###############################################################################
""" [`ArrheniusTD`](@ref) type Kc TD """
KcTDCLM(FT) =
    ArrheniusTD{FT}(  40.49, 79430.0 / GAS_R(), 79430.0 / RT_25(FT));

""" [`ArrheniusTD`](@ref) type Ko TD """
KoTDCLM(FT) =
    ArrheniusTD{FT}(27840.0, 36380.0 / GAS_R(), 36380.0 / RT_25(FT));

""" [`ArrheniusTD`](@ref) type Kpep TD """
KpepTDCLM(FT) =
    ArrheniusTD{FT}(    8.0, 36000.0 / GAS_R(), 36000.0 / RT_25(FT));

""" [`ArrheniusTD`](@ref) type Γ* TD """
ΓStarTDCLM(FT) =
    ArrheniusTD{FT}(  4.275, 37830.0 / GAS_R(), 37830.0 / RT_25(FT));

""" [`ArrheniusPeakTD`](@ref) type Jmax TD """
function JmaxTDBernacchi(FT)
    ΔHa_to_RT25::FT = 57500.0  / RT_25(FT);
    ΔHd_to_R::FT = 439000.0 / GAS_R();
    ΔSv_to_R::FT = 1400.0   / GAS_R();
    C::FT = 1 + exp( ΔSv_to_R - ΔHd_to_R/T_25(FT) );
    return ArrheniusPeakTD{FT}(ΔHa_to_RT25, ΔHd_to_R, ΔSv_to_R, C)
end

""" [`ArrheniusPeakTD`](@ref) type Jmax TD """
function JmaxTDCLM(FT)
    ΔHa_to_RT25::FT = 43540.0  / RT_25(FT);
    ΔHd_to_R::FT = 150000.0 / GAS_R();
    ΔSv_to_R::FT = 490.0    / GAS_R();
    C::FT = 1 + exp( ΔSv_to_R - ΔHd_to_R/T_25(FT) );
    return ArrheniusPeakTD{FT}(ΔHa_to_RT25, ΔHd_to_R, ΔSv_to_R, C)
end

""" [`ArrheniusPeakTD`](@ref) type Respiration TD """
function RespirationTDCLM(FT)
    ΔHa_to_RT25::FT = 46390.0  / RT_25(FT);
    ΔHd_to_R::FT = 150650.0 / GAS_R();
    ΔSv_to_R::FT = 490.0    / GAS_R();
    C::FT = 1 + exp( ΔSv_to_R - ΔHd_to_R/T_25(FT) );
    return ArrheniusPeakTD{FT}(ΔHa_to_RT25, ΔHd_to_R, ΔSv_to_R, C)
end

""" [`ArrheniusPeakTD`](@ref) type Vcmax TD """
function VcmaxTDCLM(FT)
    ΔHa_to_RT25::FT = 65330.0  / RT_25(FT);
    ΔHd_to_R::FT = 150000.0 / GAS_R();
    ΔSv_to_R::FT = 490.0    / GAS_R();
    C::FT = 1 + exp( ΔSv_to_R - ΔHd_to_R/T_25(FT) );
    return ArrheniusPeakTD{FT}(ΔHa_to_RT25, ΔHd_to_R, ΔSv_to_R, C)
end








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
FluorescenceVanDerTol(FT) =
        FluoParaSet{FT}(2.48, 2.83, 0.114);

""" [`FluoParaSet`](@ref) type parameter set using Flexas's data (drought) """
FluorescenceVanDerTolDrought(FT) =
        FluoParaSet{FT}(5.01, 1.93, 10.0);








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
