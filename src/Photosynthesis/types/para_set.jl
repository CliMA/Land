#=
The structure tree of AbstractPhotoModelParaSet
AbstractPhotoModelParaSet
---> C3ParaSet
    ---> C3ParaSetVcJ      # Vc, and J model for C3 plants
        ---> C3VcJBernacchi
    ---> C3ParaSetVcVpJ    # Vc, Vp, and J model for C3 plants
        ---> C3VcVpJBernacchi
---> C4ParaSet
    ---> C4ParaSetVcVpJ    # Vc, Vp, and J model for C4 plants
        ---> C4VcVpJBoyd
        ---> C4VcVpJCLM
=#
abstract type AbstractPhotoModelParaSet end

abstract type C3ParaSet <: AbstractPhotoModelParaSet end
abstract type C4ParaSet <: AbstractPhotoModelParaSet end

abstract type C3ParaSetVcJ   <: C3ParaSet end
abstract type C3ParaSetVcVpJ <: C3ParaSet end
abstract type C4ParaSetVcVpJ <: C4ParaSet end




"""
    C3VcJBernacchi{FT} <: C3ParaSetVcJ

A C3ParaSetVcJ type of parameter sets for different temperature dependencies.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct C3VcJBernacchi{FT} <: C3ParaSetVcJ
    "Jmax temperature dependency"
    JT ::JmaxTDBernacchi        = JmaxTDBernacchi{FT}()
    "Kc temperature dependency"
    KcT::KcTDBernacchi          = KcTDBernacchi{FT}()
    "Ko temperature dependency"
    KoT::KoTDBernacchi          = KoTDBernacchi{FT}()
    "Respiration temperature dependency"
    ReT::RespirationTDBernacchi = RespirationTDBernacchi{FT}()
    "Vcmax temperature dependency"
    VcT::VcmaxTDBernacchi       = VcmaxTDBernacchi{FT}()
    "Γ_star temperature dependency"
    ΓsT::ΓStarTDBernacchi       = ΓStarTDBernacchi{FT}()
    "Vcmax25 and respiration correlation"
    VR ::VtoRDefault            = VtoRDefault{FT}()
end




"""
    C3VcVpJBernacchi{FT} <: C3ParaSetVcVpJ

A C3ParaSetVcVpJ type of parameter sets for different temperature dependencies.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct C3VcVpJBernacchi{FT} <: C3ParaSetVcVpJ
    "Jmax temperature dependency"
    JT ::JmaxTDBernacchi        = JmaxTDBernacchi{FT}()
    "Kc temperature dependency"
    KcT::KcTDBernacchi          = KcTDBernacchi{FT}()
    "Ko temperature dependency"
    KoT::KoTDBernacchi          = KoTDBernacchi{FT}()
    "Respiration temperature dependency"
    ReT::RespirationTDBernacchi = RespirationTDBernacchi{FT}()
    "Vcmax temperature dependency"
    VcT::VcmaxTDBernacchi       = VcmaxTDBernacchi{FT}()
    "Γ_star temperature dependency"
    ΓsT::ΓStarTDBernacchi       = ΓStarTDBernacchi{FT}()
    "Vcmax25 and respiration correlation"
    VR ::VtoRDefault            = VtoRDefault{FT}()
end




"""
    C4VcVpJBoyd{FT} <: C4ParaSetVcVpJ

A C4ParaSetVcVpJ type of parameter sets for different temperature dependencies.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct C4VcVpJBoyd{FT} <: C4ParaSetVcVpJ
    "Kpep temperature dependency"
    KpT::KpepTDBoyd       = KpepTDBoyd{FT}()
    "Respiration temperature dependency"
    ReT::RespirationTDCLM = RespirationTDCLM{FT}()
    "Vcmax temperature dependency"
    VcT::VcmaxTDCLM       = VcmaxTDCLM{FT}()
    "Vpmax temperature dependency"
    VpT::VpmaxTDBoyd      = VpmaxTDBoyd{FT}()
    "Vcmax25 and respiration correlation"
    VR ::VtoRDefault      = VtoRDefault{FT}()
end




"""
    C4VcVpJCLM{FT} <: C4ParaSetVcVpJ

A C4ParaSetVcVpJ type of parameter sets for different temperature dependencies.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct C4VcVpJCLM{FT} <: C4ParaSetVcVpJ
    "Kpep temperature dependency"
    KpT::KpepTDCLM        = KpepTDCLM{FT}()
    "Respiration temperature dependency"
    ReT::RespirationTDCLM = RespirationTDCLM{FT}()
    "Vcmax temperature dependency"
    VcT::VcmaxTDCLM       = VcmaxTDCLM{FT}()
    "Vpmax temperature dependency"
    VpT::VpmaxTDBoyd      = VpmaxTDBoyd{FT}()
    "Vcmax25 and respiration correlation"
    VR ::VtoRDefault      = VtoRDefault{FT}()
end
