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
    E1  = FT(4);
    E2  = FT(8);
    return C3ParaSet{FT}(JT, KcT, KoT, ReT, VcT, ΓsT, Flu, E1, E2)
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
    E1  = FT(4);
    E2  = FT(8);
    return C3ParaSet{FT}(JT, KcT, KoT, ReT, VcT, ΓsT, Flu, E1, E2)
end

""" [`C4ParaSet`](@ref) type C4 photosynthesis using CLM5's data """
function C4CLM(FT)
    KpT = KpepTDCLM(FT);
    ReT = RespirationTDCLM(FT);
    VcT = VcmaxTDCLM(FT);
    VpT = VpmaxTDBoyd(FT);
    Flu = FluorescenceVanDerTolDrought(FT);
    return C4ParaSet{FT}(KpT, ReT, VcT, VpT, Flu)
end
