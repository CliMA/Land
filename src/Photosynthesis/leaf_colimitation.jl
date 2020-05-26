"Colimitation"
abstract type AbstractColimitation end

struct  MinColimit    <: AbstractColimitation end

Base.@kwdef struct  CurvedColimit{FT} <: AbstractColimitation 
    Θ::FT = 0.995
end


function photosynthesis_colimit!(mod::MinColimit,target, A...)
    target = min(A...);
end

function photosynthesis_colimit!(mod::CurvedColimit,target,A...)
    a = A[1]
    #@show A
    for j=2:length(A)
        lower_quadratic!(mod.Θ, -(a + A[j]), a * A[j],a)
    end
    isnan(a) ? target = minimum(A) : target =a
end
