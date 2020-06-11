###############################################################################
#
# Temperature Dedenpdency Parameter sets
# Data source: Bernacchi et al. (2001)
# Improved temperature response functions for models of Rubisco‐limited
# photosynthesis
# These parameter sets passed the FT test
# These parameter are documented in the Photosynthesis page
#
###############################################################################
""" Create Kc TD from Bernacchi's data """
KcTDBernacchi(FT)          = ArrheniusTD{FT, 41.0264925, 79430.0}()
""" Create Ko TD from Bernacchi's data """
KoTDBernacchi(FT)          = ArrheniusTD{FT,   28208.88, 36380.0}()
""" Create Respiration TD from Bernacchi's data """
RespirationTDBernacchi(FT) = ArrheniusTD{FT,        Inf, 46390.0}()
""" Create Vcmax TD from Bernacchi's data """
VcmaxTDBernacchi(FT)       = ArrheniusTD{FT,        Inf, 65330.0}()
""" Create Vomax TD from Bernacchi's data """
VomaxTDBernacchi(FT)       = ArrheniusTD{FT,        Inf, 60110.0}()
""" Create Γ* TD from Bernacchi's data """
ΓStarTDBernacchi(FT)       = ArrheniusTD{FT, 4.33164375, 37830.0}()








###############################################################################
#
# Temperature Dedenpdency Parameter sets
# Data source: Boyd et al. (2001)
# Temperature responses of C4 photosynthesis: biochemical analysis of Rubisco,
# phosphoenolpyruvate carboxylase, and carbonic anhydrase in Setaria viridis
# These parameter sets passed the FT test
# These parameter are documented in the Photosynthesis page
#
###############################################################################
""" Create Kpep TD from Boyd's data """
KpepTDBoyd(FT) = ArrheniusTD{FT, 16.0, 36300.0}()
""" Create Vpmax TD from Boyd's data """
VpmaxTDBoyd(FT) = ArrheniusPeakTD{FT, 94800.0, 73300.0, 250.0}()








###############################################################################
#
# Temperature Dedenpdency Parameter sets
# Data source: Leuning (2002)
# Temperature dependence of two parameters in a photosynthesis model
# These parameter sets passed the FT test
# These parameter are documented in the Photosynthesis page
#
###############################################################################
""" Create Jmax TD from Leuning's data """
JmaxTDLeuning(FT)  = ArrheniusPeakTD{FT, 50300.0, 152044.0, 495.0}()
""" Create Vcmax TD from Leuning's data """
VcmaxTDLeuning(FT) = ArrheniusPeakTD{FT, 73637.0, 149252.0, 486.0}()








###############################################################################
#
# Temperature Dedenpdency Parameter sets
# Data source missing
# These parameter sets passed the FT test
# These parameter are documented in the Photosynthesis page
#
###############################################################################
""" Create Kc TD """
KcTDCLM(FT)    = ArrheniusTD{FT,   40.49, 79430.0}()
""" Create Ko TD """
KoTDCLM(FT)    = ArrheniusTD{FT, 27840.0, 36380.0}()
""" Create Kpep TD """
KpepTDCLM(FT)  = ArrheniusTD{FT,     8.0, 36000.0}()
""" Create Γ* TD """
ΓStarTDCLM(FT) = ArrheniusTD{FT,   4.275, 37830.0}()
""" Create Jmax TD """
JmaxTDBernacchi(FT)  = ArrheniusPeakTD{FT, 57500.0, 439000.0, 1400.0}()
""" Create Jmax TD """
JmaxTDCLM(FT)        = ArrheniusPeakTD{FT, 43540.0, 150000.0,  490.0}()
""" Create Respiration TD """
RespirationTDCLM(FT) = ArrheniusPeakTD{FT, 46390.0, 150650.0,  490.0}()
""" Create Vcmax TD """
VcmaxTDCLM(FT)       = ArrheniusPeakTD{FT, 65330.0, 150000.0,  490.0}()








###############################################################################
#
# Vcmax to R correlation
# These parameter sets passed the FT test
# These parameter are documented in the Photosynthesis page
#
###############################################################################
""" A constant of 0.01 """
VtoRCollatz(FT) = FT(0.010)
""" A constant of 0.015 """
VtoRDefault(FT) = FT(0.015)








###############################################################################
#
# Pre-set photosynthesis model parameter sets
# These parameter sets passed the FT test
# These parameter are documented in the Photosynthesis page
#
###############################################################################
""" Bernacchi set for C3 photosynthesis """
C3Bernacchi(FT) = C3ParaSet{
                        FT,
                        JmaxTDBernacchi(FT),
                        KcTDBernacchi(FT),
                        KoTDBernacchi(FT),
                        RespirationTDBernacchi(FT),
                        VcmaxTDBernacchi(FT),
                        ΓStarTDBernacchi(FT),
                        VtoRDefault(FT),
                        4.0,
                        8.0}()
""" CLM set for C3 photosynthesis """
C3CLM(FT) = C3ParaSet{
                        FT,
                        JmaxTDCLM(FT),
                        KcTDCLM(FT),
                        KoTDCLM(FT),
                        RespirationTDCLM(FT),
                        VcmaxTDCLM(FT),
                        ΓStarTDCLM(FT),
                        VtoRDefault(FT),
                        4.0,
                        8.0}()
""" CLM set for C4 photosynthesis """
C4CLM(FT) = C4ParaSet{
                        FT,
                        KpepTDCLM(FT),
                        RespirationTDCLM(FT),
                        VcmaxTDCLM(FT),
                        VpmaxTDBoyd(FT),
                        VtoRDefault(FT)}()








###############################################################################
#
# Fluorescence model paraset
# Data source: Flexas et al. (2002)
# Energy dissipation in C3 plants under drought
# This parameter set passed the FT test
# This parameter is documented in the Photosynthesis page
#
###############################################################################
""" Create Fluorescence parameter set using Flexas's data """
FluorescenceFlexas(FT) = FluoParaSet{FT, 5.01, 1.93, 10.0}()
