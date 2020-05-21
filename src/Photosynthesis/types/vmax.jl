#=
The structure tree of AbstractVmax
AbstractVmax
---> AbstractVcmax
    ---> struct VcmaxBernacchi    # using Arrhenius correction
    ---> struct VcmaxCLM          # using Arrhenius peak correction
    ---> struct VcmaxLeuning      # using Arrhenius peak correction
---> AbstractVomax
    ---> struct VomaxBernacchi    # using Arrhenius correction
---> AbstractVpmax
    ---> struct VpmaxBoyd         # using Arrhenius peak correction
=#
abstract type AbstractVmax end

abstract type AbstractVcmax <: AbstractVmax end
abstract type AbstractVomax <: AbstractVmax end
abstract type AbstractVpmax <: AbstractVmax end




"""
    VcmaxBernacchi{FT} <: AbstractVcmax

A non-mutable AbstractVcmax type `VcmaxBernacchi` that stores information for Vcmax temperature correction.
The equation used for temperature correction is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The data source for the constants can be referred from Bernacchi et al. (2001) "Improved temperature response functions for models of Rubisco‐limited photosynthesis".

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct VcmaxBernacchi{FT} <: AbstractVcmax
    "Ratio between ΔHa and R `[K]`"
    ΔHa_to_R   ::FT = FT(65330.0) / FT(GAS_R)
    "Ratio between ΔHa and R*K-25"
    ΔHa_to_RT25::FT = FT(65330.0) / FT(GAS_R*K_25)
end




"""
    VcmaxCLM{FT} <:AbstractVcmax

A non-mutable AbstractVcmax type `VcmaxCLM` that stores information for Vcmax temperature correction.
The equation used for temperature correction is `Arrhenius equation` (see Leuning (2002) "Leuning, R. "Temperature dependence of two parameters in a photosynthesis model").
The data source for the constants can be referred from "".

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct VcmaxCLM{FT} <: AbstractVcmax
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
    VcmaxLeuning{FT} <:AbstractVcmax

A non-mutable AbstractVcmax type `VcmaxLeuning` that stores information for Vcmax temperature correction.
The equation used for temperature correction is `Arrhenius equation` (see Leuning (2002) "Leuning, R. "Temperature dependence of two parameters in a photosynthesis model").
The data source for the constants can be referred from leuning (2002).

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct VcmaxLeuning{FT} <: AbstractVcmax
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
    VomaxBernacchi{FT} <: AbstractVomax

A non-mutable AbstractVomax type `VomaxBernacchi` that stores information for Vomax temperature correction.
The equation used for temperature correction is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The data source for the constants can be referred from Bernacchi et al. (2001) "Improved temperature response functions for models of Rubisco‐limited photosynthesis".

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct VomaxBernacchi{FT} <: AbstractVomax
    "Ratio between ΔHa and R `[K]`"
    ΔHa_to_R   ::FT = FT(60110.0) / FT(GAS_R)
    "Ratio between ΔHa and R*K-25"
    ΔHa_to_RT25::FT = FT(60110.0) / FT(GAS_R*K_25)
end




"""
    VpmaxBoyd{FT} <: AbstractVpmax

A non-mutable AbstractVpmax type `VpmaxBoyd` that stores information for Vomax temperature correction.
The equation used for temperature correction is Arrhenius peak correction (see Leuning (2002) "Leuning, R. "Temperature dependence of two parameters in a photosynthesis model").
The data source for the constants can be referred from Boyd et al. (2001) "Temperature responses of C4 photosynthesis: biochemical analysis of Rubisco, phosphoenolpyruvate carboxylase, and carbonic anhydrase in Setaria viridis".

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct VpmaxBoyd{FT} <: AbstractVpmax
    "Ratio between ΔHa and R*K_25"
    ΔHa_to_RT25::FT = FT(94800.0) / FT(GAS_R*K_25)
    "Ratio between ΔHd and R"
    ΔHd_to_R   ::FT = FT(73300.0) / FT(GAS_R)
    "Ratio between ΔSv and R"
    ΔSv_to_R   ::FT = FT(250.0) / FT(GAS_R)
    "Correction factor C = 1 + exp( Sv/R + Hd/(RT0) )"
    C          ::FT = 1 + exp( ΔSv_to_R - ΔHd_to_R/FT(K_25) )
end
