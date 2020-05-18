"Electron Transport Rate parameters"
abstract type AbstractJmax end
abstract type AbstractVmax end

"Carboxylation Rate Parameters"
abstract type AbstractVcmax <: AbstractVmax end

"PEP Carboxylase Rate Parameters"
abstract type AbstractVpmax <: AbstractVmax end

"Michaelis Menten (MM) Parameters"
abstract type AbstractMM end



"Jmax Rate Parameters as in CLM5"
Base.@kwdef struct JmaxCLM{FT} <: AbstractJmax
    "Activation energy for ETR (J/mol)"
    ΔHa::FT    = 43540.;                      
    "Deactivation energy for ETR (J/mol)"
    ΔHd::FT    = 150000.;
    "Entropy term for ETR (J/mol/K)"
    ΔS::FT    = 490.;                       
    "Scaling factor to ensure Jmax(25)==Jmax25 "
    scale::FT    = fth25(ΔHd, ΔS);
end

"Carboxylation Rate Parameters as in CLM5"
Base.@kwdef struct VcmaxCLM{FT} <: AbstractVcmax
    "Activation energy for ETR (J/mol)"
    ΔHa::FT    = 65330.;                      
    "Deactivation energy for ETR (J/mol)"
    ΔHd::FT    = 150000.;
    "Entropy term for ETR (J/mol/K)" # Can be adjusted with growth temperature later!
    ΔS::FT    = 490.;                       
    "Scaling factor to ensure Vcmax(25)==Vcmax25 "
    scale::FT    = fth25(ΔHd, ΔS);
end

"Jmax Rate Parameters as in CLM5"
Base.@kwdef struct JmaxBernacchi{FT} <: AbstractJmax
    "Activation energy for ETR (J/mol)"
    ΔHa::FT    = 57500.;                      
    "Deactivation energy for ETR (J/mol)"
    ΔHd::FT    = 439000.;
    "Entropy term for ETR (J/mol/K)"
    ΔS::FT    = 1400.;                       
    "Scaling factor to ensure Jmax(25)==Jmax25 "
    scale::FT    = fth25(ΔHd, ΔS);
end

"Carboxylation Rate Parameters as in CLM5"
Base.@kwdef struct VcmaxBernacchi{FT} <: AbstractVcmax
    "Activation energy for ETR (J/mol)"
    ΔHa::FT    = 65330.;                                 
    "Scaling factor to ensure Vcmax(25)==Vcmax25 "
    scale::FT    = fArrhenius(FT(298.15), ΔHa);
end

"
PEP Carboxylation Rate Parameters 
from Boyd et al, 2015, Plant Physiology: Temperature Responses of C4 Photosynthesis:
Biochemical Analysis of Rubisco, Phosphoenolpyruvate Carboxylase, and Carbonic Anhydrase in
Setaria viridis
"
Base.@kwdef struct Vpmax{FT} <: AbstractVpmax
    "Activation energy for ETR (J/mol)"
    ΔHa::FT    = 95000.;                      
    "Deactivation energy for ETR (J/mol)"
    ΔHd::FT    = 73000.;
    "Entropy term for ETR (J/mol/K)"
    ΔS::FT    = 250.;                       
    "Scaling factor to ensure Vcmax(25)==Vcmax25 "
    scale::FT    = fth25(ΔHd, ΔS);
end


"Michaelis Rate Parameters as in CLM5"
Base.@kwdef struct MM_CLM{FT} <: AbstractMM
    "Michaelis-Menten constant for CO2 at 25C (Pa)"
    Kc_25::FT = 40.49;     
    "Michaelis-Menten constant for O2 at 25C  (Pa)"
    Ko_25::FT = 27.84;
    "Michaelis-Menten constant for PEP Carboxylase at 25C  (Pa) (from von Caemmerer Book 2000) "
    Kpep_25::FT = 8.0;                   
    "Standard CO2 compensation point at 25C (Pa)"
    Γ_25::FT = 4.275;# (Can add O2 dependence later) *o₂/0.209;     
    "Activation energy for CO2 MM constant (J/mol)"
    ΔHa_Kc::FT   = 79430.;                      
    "Activation energy for O2 MM constant (J/mol)"
    ΔHa_Ko::FT   = 36380.;
    "Activation energy for CO2 compensation point (J/mol)"
    ΔHa_Γ::FT    = 37830.;
    "Activation energy for PEP MM constant (J/mol)"
    ΔHa_Kp::FT   = 36000.;                      
end


# Functions:
""" 
    max_carboxylation_rate(model::VcmaxCLM, leaf)
Calculates the maximum carboxylation rate (Vcmax) at the leaf temperature
"""
function max_carboxylation_rate!(model::VcmaxCLM, l)
    @unpack ΔHa, ΔHd, ΔS, scale = model
    l.Vcmax  = l.Vcmax25 * ft(l.T, ΔHa) * fth(l.T, ΔHd, ΔS, scale);
end

""" 
    max_carboxylation_rate(model::VcmaxBernacchi, leaf)
Calculates the maximum carboxylation rate (Vcmax) at the leaf temperature (Bernacchi 2001 paper)
"""
function max_carboxylation_rate!(model::VcmaxBernacchi, l)
    @unpack ΔHa,scale = model
    l.Vcmax  = l.Vcmax25 * fArrhenius(l.T, ΔHa) / scale;
end


""" 
    max_electron_transport_rate(model::JmaxCLM, leaf)
Calculates the potential max_electron transport rate (Jmax) at the leaf temperature 
"""
function max_electron_transport_rate!(model::AbstractJmax, l)
    @unpack ΔHa, ΔHd, ΔS, scale = model
    l.Jmax  = l.Jmax25 * ft(l.T, ΔHa) * fth(l.T, ΔHd, ΔS, scale);
end

""" 
    michaelis_menten_constants(model::JmaxCLM, leaf)
Calculates the leaf temperature adjusted Michaelis Menten constants  
"""
function michaelis_menten_constants!(model::AbstractMM, l)
    @unpack ΔHa_Kc, ΔHa_Ko, ΔHa_Kp, ΔHa_Γ, Kc_25, Ko_25, Kpep_25,Γ_25  = model
    l.Kc  = Kc_25 * ft(l.T, ΔHa_Kc);
    l.Ko  = Ko_25 * ft(l.T, ΔHa_Ko);
    l.Kpep  = Kpep_25 * ft(l.T, ΔHa_Kp);
    l.Γstar   = Γ_25 * ft(l.T, ΔHa_Γ);
end

# Scaling functions for Photosynthesis temperature response and inhibition
function fArrhenius(tl, Ea)
    FT = eltype(tl)
    Rgas = FT(8.31446261815324);
    exp(-Ea/(Rgas*tl));
end;

function ft(tl, ha)
    FT = eltype(tl)
    Rgas = FT(8.31446261815324);
    T25 = FT(298.15) 
    exp(ha/(Rgas*(T25)) * (1-(T25)/tl));
end;

function fth(tl, hd, se, fc)
    FT = eltype(tl)
    Rgas = FT(8.31446261815324);
    fc / (1 + exp((-hd+se*tl)/(Rgas*tl)));
end;

function fth25(hd, se) 
    FT = eltype(hd)
    T25 = FT(298.15);
    Rgas = FT(8.31446261815324); 
    1 + exp( (-hd + se * T25) / (Rgas * T25 ));
end;
