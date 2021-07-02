###############################################################################
#
# Photosynthesis model parameter set -- temperature dependencies and etc
#
###############################################################################
#= AbstractPhotoModelParaSet type tree
---> C3ParaSet
---> C4ParaSet
=#
"""
    abstract type AbstractPhotoModelParaSet{FT}

Hierarchy of the `AbstractPhotoModelParaSet`:
- [`C3ParaSet`](@ref)
- [`C4ParaSet`](@ref)
"""
abstract type AbstractPhotoModelParaSet{FT} end




"""
    mutable struct C3Cytochrome{FT}

Parameter sets for C3 photosynthesis with Cytochrome activity.

# Fields
$(TYPEDFIELDS)
"""
mutable struct C3Cytochrome{FT<:AbstractFloat} <: AbstractPhotoModelParaSet{FT}
    "Coefficient 4.0/4.5 for NADPH/ATP requirement stochiometry, respectively"
    Eff_1::FT
    "Coefficient 8.0/10.5 for NADPH/ATP requirement stochiometry, respectively"
    Eff_2::FT
end




"""
    mutable struct C3Paraset{FT}

Parameter sets for C3 photosynthesis.

# Fields
$(TYPEDFIELDS)
"""
mutable struct C3ParaSet{FT<:AbstractFloat} <: AbstractPhotoModelParaSet{FT}
    "Jmax temperature dependency"
    JT::AbstractTDParameterSet{FT}
    "Kc temperature dependency"
    KcT::AbstractTDParameterSet{FT}
    "Ko temperature dependency"
    KoT::AbstractTDParameterSet{FT}
    "Respiration temperature dependency"
    ReT::AbstractTDParameterSet{FT}
    "Vcmax temperature dependency"
    VcT::AbstractTDParameterSet{FT}
    "Γ_star temperature dependency"
    ΓsT::AbstractTDParameterSet{FT}
    "Fluorescence model"
    Flu::AbstractFluoModelParaSet{FT}
    "Vcmax25 and respiration correlation"
    VR::FT
    "Coefficient 4.0/4.5 for NADPH/ATP requirement stochiometry, respectively"
    Eff_1::FT
    "Coefficient 8.0/10.5 for NADPH/ATP requirement stochiometry, respectively"
    Eff_2::FT
end




"""
    mutable struct C4ParaSet{FT}

Parameter sets for C3 photosynthesis.

# Fields
$(TYPEDFIELDS)
"""
mutable struct C4ParaSet{FT<:AbstractFloat} <: AbstractPhotoModelParaSet{FT}
    "Kpep temperature dependency"
    KpT::AbstractTDParameterSet{FT}
    "Respiration temperature dependency"
    ReT::AbstractTDParameterSet{FT}
    "Vcmax temperature dependency"
    VcT::AbstractTDParameterSet{FT}
    "Vpmax temperature dependency"
    VpT::AbstractTDParameterSet{FT}
    "Fluorescence model"
    Flu::AbstractFluoModelParaSet{FT}
    "Vcmax25 and respiration correlation"
    VR::FT
end
