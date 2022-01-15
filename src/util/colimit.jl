
abstract type AbstractColimit{FT<:AbstractFloat} end


struct MinimumColimit{FT<:AbstractFloat} <:AbstractColimit{FT} end


mutable struct QuadraticColimit{FT<:AbstractFloat} <: AbstractColimit{FT}
    CURVATURE::FT
end
