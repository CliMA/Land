abstract type AbstractLeafRespiration end

"""
"""
Base.@kwdef struct RespirationCLM{FT} <: AbstractLeafRespiration
    "Activation energy for Respiration `[J/mol]`"
    Ea::FT    = 46390.;                      
    "Deactivation energy for Respiration `[J/mol]`"
    Hd::FT    = 150650.;
    "Entropy term for Respiration `[J/mol/K]`"
    ΔS::FT    = 490.;                       # Entropy term for Rd (J/mol/K)
    "Scaling factor for high temperature inhibition (25 C = 1.0)"
    scale::FT    = fth25(Hd, ΔS);          # Scaling factor for high temperature inhibition (25 C = 1.0)
end

Base.@kwdef struct RespirationBernacchi{FT} <: AbstractLeafRespiration
    "Activation energy for Respiration (J/mol)"
    Ea::FT    = 46390.;                      
end


function leaf_respiration!(model::RespirationCLM, l)
    @unpack Ea, Hd, ΔS, scale = model
    l.Rd  = l.Rd25 * ft(l.T, Ea) * fth(l.T, Hd, ΔS, scale);
end

function leaf_respiration!(model::RespirationBernacchi, l)
    @unpack Ea = model
    l.Rd  = l.Rd25 * ft(l.T, Ea);
end