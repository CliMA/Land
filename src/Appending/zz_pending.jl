# What is this?

"Colimitation"
abstract type AbstractColimitation end

struct  MinColimit    <: AbstractColimitation end

Base.@kwdef struct  CurvedColimit{FT} <: AbstractColimitation 
    Θ::FT = 0.995
end


function photosynthesis_colimit!(mod::MinColimit,A...)
    min(A...);
end

function photosynthesis_colimit!(mod::CurvedColimit,A...)
    a = A[1]
    #@show A
    for j=2:length(A)
        a = lower_quadratic(mod.Θ, -(a + A[j]), a * A[j])
    end
    #@show a
    isnan(a) ? minimum(A) : a
    #@show target
end





#=
"""
PhotoMods
describes all necessary modules
"""
Base.@kwdef struct PhotoMods{FM,PM,RM,SM,JM,VM,MM,BL,CL} <: AbstractPhotosynthesis
    fluorescence::FM     = FlexasTolBerryFluorescence{FT}()
    photosynthesis::PM   = C3FvCBPhoto()
    respiration::RM      = RespirationCLM{FT}()
    stomatal::SM         = BallBerryStomata{FT}()
    Jmax::JM             = JmaxCLM{FT}()
    Vmax::VM             = VcmaxCLM{FT}()
    MichaelisMenten::MM  = MM_CLM{FT}()
    BoundaryLayer::BL    = FixedBoundaryResistance{FT}(ra=1)
    colimitation::CL     = CurvedColimit{FT}()
end
=#
