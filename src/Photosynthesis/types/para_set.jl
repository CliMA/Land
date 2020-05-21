#=
The structure tree of AbstractPhotoModelParaSet
AbstractPhotoModelParaSet
---> PhotoModelParaSetVcJ      # Vc, and J model
---> PhotoModelParaSetVcVpJ    # Vc, Vp, and J model
=#
abstract type AbstractPhotoModelParaSet end

abstract type PhotoModelParaSetVcJ   <: AbstractPhotoModelParaSet end
abstract type PhotoModelParaSetVcVpJ <: AbstractPhotoModelParaSet end




"""
    PMVcJBernacchi{FT} <: PhotoModelParaSetVcJ

A PhotoModelParaSetVcJ type of parameter sets for different temperature dependencies.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct PMPSVcJBernacchi{FT} <: PhotoModelParaSetVcJ
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
end
