#=
The structure tree of AbstractVmaxTD
AbstractVmaxTD
---> AbstractVcmaxTD
    ---> struct VcmaxTDBernacchi    # using Arrhenius correction
    ---> struct VcmaxTDCLM          # using Arrhenius peak correction
    ---> struct VcmaxTDLeuning      # using Arrhenius peak correction
---> AbstractVomaxTD
    ---> struct VomaxBernacchi    # using Arrhenius correction
---> AbstractVpmaxTD
    ---> struct VpmaxBoyd         # using Arrhenius peak correction
=#
abstract type AbstractVmaxTD end

abstract type AbstractVcmaxTD <: AbstractVmaxTD end
abstract type AbstractVomaxTD <: AbstractVmaxTD end
abstract type AbstractVpmaxTD <: AbstractVmaxTD end




"""
    VcmaxTDBernacchi{FT} <: AbstractVcmaxTD

A non-mutable AbstractVcmaxTD type `VcmaxTDBernacchi` that stores information for Vcmax temperature correction.
The equation used for temperature correction is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The data source for the constants can be referred from Bernacchi et al. (2001) "Improved temperature response functions for models of Rubisco‐limited photosynthesis".

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct VcmaxTDBernacchi{FT} <: AbstractVcmaxTD
    "Ratio between ΔHa and R `[K]`"
    ΔHa_to_R   ::FT = FT(65330.0) / FT(GAS_R)
    "Ratio between ΔHa and R*K-25"
    ΔHa_to_RT25::FT = FT(65330.0) / FT(GAS_R*K_25)
end




"""
    VcmaxTDCLM{FT} <:AbstractVcmaxTD

A non-mutable AbstractVcmaxTD type `VcmaxTDCLM` that stores information for Vcmax temperature correction.
The equation used for temperature correction is `Arrhenius equation` (see Leuning (2002) "Leuning, R. "Temperature dependence of two parameters in a photosynthesis model").
The data source for the constants can be referred from "".

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct VcmaxTDCLM{FT} <: AbstractVcmaxTD
    "Ratio between ΔHa and R*K_25"
    ΔHa_to_RT25::FT = FT(65330.0) / FT(GAS_R*K_25)
    "Ratio between ΔHd and R"
    ΔHd_to_R   ::FT = FT(150000.0) / FT(GAS_R)
    "Ratio between ΔSv and R"
    ΔSv_to_R   ::FT = FT(490.0) / FT(GAS_R)
    "Correction factor C = 1 + exp( Sv/R + Hd/(RT0) )"
    C          ::FT = 1 + exp( ΔSv_to_R - ΔHd_to_R/FT(K_25) )
end




"""
    VcmaxTDLeuning{FT} <:AbstractVcmaxTD

A non-mutable AbstractVcmaxTD type `VcmaxTDLeuning` that stores information for Vcmax temperature correction.
The equation used for temperature correction is `Arrhenius equation` (see Leuning (2002) "Leuning, R. "Temperature dependence of two parameters in a photosynthesis model").
The data source for the constants can be referred from leuning (2002).

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct VcmaxTDLeuning{FT} <: AbstractVcmaxTD
    "Ratio between ΔHa and R*K_25"
    ΔHa_to_RT25::FT = FT(73637.0) / FT(GAS_R*K_25)
    "Ratio between ΔHd and R"
    ΔHd_to_R   ::FT = FT(149252.0) / FT(GAS_R)
    "Ratio between ΔSv and R"
    ΔSv_to_R   ::FT = FT(486.0) / FT(GAS_R)
    "Correction factor C = 1 + exp( Sv/R + Hd/(RT0) )"
    C          ::FT = 1 + exp( ΔSv_to_R - ΔHd_to_R/FT(K_25) )
end




"""
    VomaxTDBernacchi{FT} <: AbstractVomaxTD

A non-mutable AbstractVomaxTD type `VomaxTDBernacchi` that stores information for Vomax temperature correction.
The equation used for temperature correction is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The data source for the constants can be referred from Bernacchi et al. (2001) "Improved temperature response functions for models of Rubisco‐limited photosynthesis".

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct VomaxTDBernacchi{FT} <: AbstractVomaxTD
    "Ratio between ΔHa and R `[K]`"
    ΔHa_to_R   ::FT = FT(60110.0) / FT(GAS_R)
    "Ratio between ΔHa and R*K-25"
    ΔHa_to_RT25::FT = FT(60110.0) / FT(GAS_R*K_25)
end




"""
    VpmaxTDBoyd{FT} <: AbstractVpmaxTD

A non-mutable AbstractVpmaxTD type `VpmaxTDBoyd` that stores information for Vomax temperature correction.
The equation used for temperature correction is Arrhenius peak correction (see Leuning (2002) "Leuning, R. "Temperature dependence of two parameters in a photosynthesis model").
The data source for the constants can be referred from Boyd et al. (2001) "Temperature responses of C4 photosynthesis: biochemical analysis of Rubisco, phosphoenolpyruvate carboxylase, and carbonic anhydrase in Setaria viridis".

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct VpmaxTDBoyd{FT} <: AbstractVpmaxTD
    "Ratio between ΔHa and R*K_25"
    ΔHa_to_RT25::FT = FT(94800.0) / FT(GAS_R*K_25)
    "Ratio between ΔHd and R"
    ΔHd_to_R   ::FT = FT(73300.0) / FT(GAS_R)
    "Ratio between ΔSv and R"
    ΔSv_to_R   ::FT = FT(250.0) / FT(GAS_R)
    "Correction factor C = 1 + exp( Sv/R + Hd/(RT0) )"
    C          ::FT = 1 + exp( ΔSv_to_R - ΔHd_to_R/FT(K_25) )
end
