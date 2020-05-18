abstract type AbstractStomatalModel end

# Empirical Classical Models:
abstract type AbstractClassicalStomatalModel <: AbstractStomatalModel end

# Ball Berry or Medlyn model
#    g0::TT = 0.1;                              # Ball-Berry minimum leaf conductance, unstressed (mol/m2/s) - Sellers  et  al.  (1996)  used  0.1 b=  for  C3  plants  and  0.4 b=for  C4  plants
#    g1_BB::TT = 9.0;                           # Ball-Berry slope of conductance-photosynthesis relationship, unstressed - m=  for  C3  plants  and  4m=  for  C4  plants  (Collatz  et  al.  1991,  1992)
#    g1_Medlyn::TT = 126.49;                    # Medlyn slope of conductance-photosynthesis relationship, unstressed - Pa^(1/2) Medlyn et al. 2017


Base.@kwdef struct BallBerryStomata{FT} <: AbstractClassicalStomatalModel
    "Ball-Berry minimum leaf conductance (moles/m2/s)"
    g0::FT = 0.025;                      
    "Ball-Berry slope of conductance-photosynthesis relationship"
    g1::FT = 9.0;
end

Base.@kwdef struct GentineStomata{FT} <: AbstractClassicalStomatalModel
    "Ball-Berry minimum leaf conductance (moles/m2/s)"
    g0::FT = 0.025;                      
    "Ball-Berry slope of conductance-photosynthesis relationship"
    g1::FT = 9.0;
end

Base.@kwdef struct MedlynStomata{FT} <: AbstractLeafRespiration
    "Medlyn slope of conductance-photosynthesis relationship, unstressed - Pa^(1/2) Medlyn et al. 2017"
    g1::FT = 125;
    "Minimum leaf conductance,"
    g0::FT = 0.1;                        
end

" Ball-Berry stomatal conductance model:"
function stomatal_conductance!(mod::BallBerryStomata, l)
  #  Cs  : CO2 at leaf surface [ppm]
  #  RH  : relative humidity [0-1]
  #  An   : Net assimilation in 'same units of CO2 as Cs' micromoles/m2/s
  #  gs   : moles/m2/s
  gs = mod.g1 * max(l.An,1e-9) * l.RH/l.Cs  + mod.g0;
  l.dynamic_state ?   l.gs_ss = gs : l.gs = gs
end # function


" Medlyn stomatal conductance model:"
function stomatal_conductance!(mod::MedlynStomata, l)
  #  Cs  : CO2 at leaf surface
  #  VPD  : vapor pressure deficit - Pa
  #  Cs  : CO2 at leaf surface [ppm]
  #  RH  : relative humidity [0-1]
  #  An   : Net assimilation in 'same units of CO2 as Cs' micromoles/m2/s
  #  gs   : moles/m2/s
  gs = (1 +mod.g1/sqrt(l.VPD)) * max(l.An,1e-9) /l.Cs  + mod.g0;
  l.dynamic_state ?   l.gs_ss = gs : l.gs = gs 
end # function


"Gentine stomatal conductance model:"
function stomatal_conductance!(mod::GentineStomata, l)
  #  Cs  : CO2 at leaf surface [ppm]
  #  RH  : relative humidity [0-1]
  #  An   : Net assimilation in 'same units of CO2 as Cs' micromoles/m2/s
  #  gs   : moles/m2/s
  setLeafkl!(l, l.psi_l) # set hydraulic conductivity of leaf
  gs = mod.g1*l.kleaf/l.kmax * max(l.An,1e-9) /l.Cs  + mod.g0;
  l.dynamic_state ?   l.gs_ss = gs : l.gs = gs
end # function

