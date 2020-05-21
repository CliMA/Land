#=
The structure tree of AbstractJmaxTD
AbstractJmaxTD
---> struct JmaxTDBernacchi    # using Arrhenius peak correction
---> struct JmaxTDCLM          # using Arrhenius peak correction
---> struct JmaxTDLeuning      # using Arrhenius peak correction
=#
abstract type AbstractJmaxTD end




"""
    JmaxTDBernacchi{FT} <:AbstractJmaxTD

A non-mutable AbstractJmaxTD type `JmaxTDBernacchi` that stores information for Jmax temperature correction.
The equation used for temperature correction is `Arrhenius equation` (See Leuning (2002) "Leuning, R. "Temperature dependence of two parameters in a photosynthesis model").
The data source for the constants can be referred from "".

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct JmaxTDBernacchi{FT} <: AbstractJmaxTD
    "Ratio between ΔHa and R*K_25"
    ΔHa_to_RT25::FT = FT(57500.0) / FT(GAS_R*K_25)
    "Ratio between ΔHd and R"
    ΔHd_to_R   ::FT = FT(439000.0) / FT(GAS_R)
    "Ratio between ΔSv and R"
    ΔSv_to_R   ::FT = FT(1400.0) / FT(GAS_R)
    "Correction factor C = 1 + exp( Sv/R + Hd/(RT0) )"
    C          ::FT = 1 + exp( ΔSv_to_R - ΔHd_to_R/FT(K_25) )
end




"""
    JmaxTDCLM{FT} <:AbstractJmaxTD

A non-mutable AbstractJmaxTD type `JmaxTDCLM` that stores information for Jmax temperature correction.
The equation used for temperature correction is `Arrhenius equation` (See Leuning (2002) "Leuning, R. "Temperature dependence of two parameters in a photosynthesis model").
The data source for the constants can be referred from "".

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct JmaxTDCLM{FT} <: AbstractJmaxTD
    "Ratio between ΔHa and R*K_25"
    ΔHa_to_RT25::FT = FT(43540.0) / FT(GAS_R*K_25)
    "Ratio between ΔHd and R"
    ΔHd_to_R   ::FT = FT(150000.0) / FT(GAS_R)
    "Ratio between ΔSv and R"
    ΔSv_to_R   ::FT = FT(490.0) / FT(GAS_R)
    "Correction factor C = 1 + exp( Sv/R + Hd/(RT0) )"
    C          ::FT = 1 + exp( ΔSv_to_R - ΔHd_to_R/FT(K_25) )
end




"""
    JmaxTDLeuning{FT} <:AbstractJmaxTD

A non-mutable AbstractJmaxTD type `JmaxTDLeuning` that stores information for Jmax temperature correction.
The equation used for temperature correction is `Arrhenius equation` (See Leuning (2002) "Leuning, R. "Temperature dependence of two parameters in a photosynthesis model").
The data source for the constants can be referred from Leuning (2002).

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct JmaxTDLeuning{FT} <: AbstractJmaxTD
    "Ratio between ΔHa and R*K_25"
    ΔHa_to_RT25::FT = FT(50300.0) / FT(GAS_R*K_25)
    "Ratio between ΔHd and R"
    ΔHd_to_R   ::FT = FT(152044.0) / FT(GAS_R)
    "Ratio between ΔSv and R"
    ΔSv_to_R   ::FT = FT(495.0) / FT(GAS_R)
    "Correction factor C = 1 + exp( Sv/R + Hd/(RT0) )"
    C          ::FT = 1 + exp( ΔSv_to_R - ΔHd_to_R/FT(K_25) )
end
