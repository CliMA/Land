###############################################################################
#
# Temperature dependency parameter set
#
###############################################################################
#= AbstractTDParameterSet type tree
---> ArrheniusTD
---> ArrheniusPeakTD
=#
"""
    abstract type AbstractTDParameterSet{FT}

Hierarchy of the `AbstractTDParameterSet`:
- [`ArrheniusTD`](@ref)
- [`ArrheniusPeakTD`](@ref)
"""
abstract type AbstractTDParameterSet{FT} end




"""
    struct ArrheniusTD{FT}

An [`AbstractTDParameterSet`](@ref) type struct using
```math
corr = \\exp \\left( \\dfrac{ΔHa}{R T_0} - \\dfrac{ΔHa}{R T_1} \\right)
```

# Fields
$(TYPEDFIELDS)
"""
struct ArrheniusTD{FT<:AbstractFloat} <: AbstractTDParameterSet{FT}
    "Uncorrected value at 298.15 K"
    VAL_25::FT
    "Ratio between ΔHa and R `[K]`"
    ΔHa_to_R::FT
    "Ratio between ΔHa and R*K_25"
    ΔHa_to_RT25::FT
end




"""
    struct ArrheniusPeakTD{FT}

An [`AbstractTDParameterSet`](@ref) type struct using
```math
corr = \\exp \\left( \\dfrac{ΔHa}{R T_0} - \\dfrac{ΔHa}{R T_1} \\right)
       \\cdot
       \\dfrac{ 1 + \\exp \\left( \\dfrac{S_v T_0 - H_d}{R T_0} \\right) }
              { 1 + \\exp \\left( \\dfrac{S_v T_1 - H_d}{R T_1} \\right) }
```

# Fields
$(TYPEDFIELDS)
"""
struct ArrheniusPeakTD{FT<:AbstractFloat} <: AbstractTDParameterSet{FT}
    "Ratio between ΔHa and R*K_25"
    ΔHa_to_RT25::FT
    "Ratio between ΔHd and R"
    ΔHd_to_R::FT
    "Ratio between ΔSv and R"
    ΔSv_to_R::FT
    "Correction factor C = 1 + exp( Sv/R + Hd/(RT0) )"
    C::FT
end




"""
    struct Q10TD{FT}

An [`AbstractTDParameterSet`](@ref) type struct using
```math
VAL = VAL_{REF} \\left( \\dfrac{T_1 - T_{REF}}}{10} \\right)^{Q_{10}}
```

# Fields
$(TYPEDFIELDS)
"""
struct Q10TD{FT<:AbstractFloat} <: AbstractTDParameterSet{FT}
    "Uncorrected value at reference temperature"
    VAL_REF::FT
    "Reference temperature `[K]`"
    T_REF::FT
    "Power of Q10 correction"
    Q_10::FT
end
