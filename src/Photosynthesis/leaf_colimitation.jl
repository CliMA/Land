"Colimitation"
abstract type AbstractColimitation end

struct  MinColimit    <: AbstractColimitation end

Base.@kwdef struct  CurvedColimit{FT} <: AbstractColimitation 
    Θ::FT = 0.995
end


function photosynthesis_colimit(mod::MinColimit, A...)
    min(A...)
end

function photosynthesis_colimit(mod::CurvedColimit, A...)
    a = A[1]
    #@show A
    for j=2:length(A)
        a = lower_quadratic(mod.Θ, -(a + A[j]), a * A[j])
    end
    isnan(a) ? minimum(A) : a
end
