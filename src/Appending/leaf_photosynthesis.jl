abstract type AbstractPhotosynthesisModel end

abstract type AbstractC3Photosynthesis <: AbstractPhotosynthesisModel end
"""
Classical C3 Photosynthesis model from FvCB, but using a pre-defined gs value as constraint
"""
struct  C3FvCBPhotoGs  <: AbstractC3Photosynthesis end





############  model::C3FvCBPhotoGs #################################
"""
    rubisco_limited_rate!(model::C3FvCBPhotoGs, leaf, met)
Solution when Rubisco activity is limiting and a given gs
Solves quadratic equation 
"""
function rubisco_limited_rate!(model::C3FvCBPhotoGs, leaf, met)
    @unpack  Γstar, Cc, Kc, Ko, o₂, Vcmax, gs, gm, Rd, Γstar = leaf
    @unpack  Ca,ra,g_m_s_to_mol_m2_s = met
    d = Kc*(1+met.ppm_to_Pa*1e6*o₂/Ko)
    a = (ra/g_m_s_to_mol_m2_s + 1.6/gs + 1.0/gm) # = 1/gleaf
    b = -(Ca + d) - (Vcmax-Rd)*a
    c = Vcmax * (Ca - Γstar) - Rd * (Ca+d)
    leaf.Ac = lower_quadratic(a,b,c)
end

"""
    light_limited_rate!(model::C3FvCBPhotoGs,  leaf, met, APAR)
Solution when Electron Transport Rate is limiting and a given gs
C3: RuBP-limited photosynthesis (this is the NADPH requirement stochiometry)
Using curvature θ_j and Jmax
"""
function light_limited_rate!(model::C3FvCBPhotoGs,  leaf, met, APAR)
    @unpack  Γstar, Cc, gs, gm, Rd = leaf
    @unpack  Ca, ra, g_m_s_to_mol_m2_s = met
    electron_transport_rate!(model, leaf, APAR)
    # Check out Bonan's book, this equates diffusion limited and demand limited 
    a = 4*(ra/g_m_s_to_mol_m2_s + 1.6/gs + 1.0/gm) # = 1/gleaf
    b = -(4Ca + 8Γstar) - (leaf.Je - 4Rd)*a/4
    c = leaf.Je * (Ca-Γstar) - Rd * (4Ca + 8Γstar)
    leaf.Aj = lower_quadratic(a,b,c)
end









# go to plant hydraulics
# Need to ba abstracted as well!
function setkx!(l::LeafPhotosynthesisParams, psis, psi_l) # set hydraulic conductivitytimes Delta Psi
    l.kx = l.kmax * IntWeibull(psis,psi_l,l.psi_l50,l.ck)/max(psis-psi_l,1e-6); # kmax . int_psis^psil k(x)dx = kmax . IntWeibull(psil);
    #println("k_xylem = ",l.kx," psi_s=",psis," psi_l=",psi_l)
end

function setLeafkl!(l::LeafPhotosynthesisParams, psi_l) # set hydraulic conductivity
    l.kleaf = l.kmax * Weibull(psi_l,l.psi_l50,l.ck); # kmax . int_psis^psil k(x)dx = kmax . IntWeibull(psil);
end














# Go to SPAC model

####
#### Surface energy balance
####


abstract type AbstractPhotosynthesis end

export LeafEnergyWaterBalance






"""
    LeafPhotosynthesis!(mod::AbstractPhotosynthesis, leaf::leaf_params, met::meteo)

Compute net assimilation rate A, fluorescence F using biochemical model, given
- `mod` One [`PhotoMods`](@ref) model setup structure
- `leaf` One [`leaf_params`](@ref) structure 
- `met` One [`lmeteo`](@ref) structure 

"""
function LeafPhotosynthesis!(mo::AbstractPhotosynthesis, leaf::leaf_params, met::meteo)
  # Set leaf temperature:
  set_leaf_temperature!(mo, leaf)
  met.g_m_s_to_mol_m2_s = (met.P_air-met.e_air)/(physcon.Rgas*met.T_air) 
  met.ppm_to_Pa = (met.P_air-met.e_air)*1e-6
  isnan(leaf.gleaf) ?  leaf.gleaf = 0.01 :
  isnan(leaf.gs) ?  leaf.gs = 0.01 :

  # Compute Cc and Photosynthesis:
  if leaf.dynamic_state
    CcFunc!(leaf.Cc; mods=mo, leaf=leaf, met=met)
  else
    #@show (CcFunc!(leaf.Cc; mods=mo, leaf=leaf, met=met))
    @inline f(x) = CcFunc!(x; mods=mo, leaf=leaf, met=met)
    sol = find_zero(f, SecantMethod{FT}(met.Ca*0.6, met.Ca), CompactSolution(),SolutionTolerance{FT}(1e-3))
    leaf.Cc = sol.root
    #@show sol
    #@show leaf.Cc
    #CcFunc!(mo, leaf, met)
    #@show leaf.Cc
  end

  # Compute Fluorescence
  leaf_fluorescence!(mo.fluorescence,  leaf);
  return leaf.An

end # LeafPhotosynthesis (similar to biochem in SCOPE)


"""
    CcFunc!(mods::PhotoMods,  leaf::leaf_params, met::meteo)

Compute gross and net assimilation rates Ag, An using biochemical model and given Cc, using
- `model` One [`PhotoMods`](@ref) model setup structure
- `leaf` One [`leaf_params`](@ref) structure 
- `met` One [`lmeteo`](@ref) structure 

"""
function CcFunc!(Cc; mods::AbstractPhotosynthesis,  leaf::leaf_params, met::meteo)
    leaf.Cc = max(0,Cc)
    #@show leaf.Cc
    # Compute Rubisco Limited Photosynthesis
    rubisco_limited_rate!(mods.photosynthesis, leaf, met)
    # Light limited rate
    light_limited_rate!(mods.photosynthesis, leaf, met, leaf.APAR)
    # Product limited rate
    product_limited_rate!(mods.photosynthesis, leaf, met)
    
    # adjust aerodynamic resistance based on leaf boundary layer and Monin Obukhov
    boundary_layer_resistance!(mods.BoundaryLayer, leaf,  met)

    @unpack Ac, Aj,Ap, Cs, VPD, RH, gm, gs, Ag, Rd, An, gleaf, esat, Cc = leaf
    #@show Ac, Aj,Ap, Cs, VPD, RH, gm, gs, Ag, Rd, An, gleaf, esat, Cc 
    @unpack Ca, ra, g_m_s_to_mol_m2_s, e_air = met

    # add colimitation later:
    #@show Ac,Aj,Ap, Cc
    leaf.Ag = photosynthesis_colimit!(mods.colimitation, Ac, Aj, Ap) 

    # Net photosynthesis due to biochemistry
    leaf.An = leaf.Ag - leaf.Rd # net assimilation
    #@show leaf.An, leaf.Ag
    # Compute stomatal conductance (update gs directly or is gs_ss if model is dynamic!):
    stomatal_conductance!(mods.stomatal, leaf)

    
    #@show leaf.Ag
    # CO2 at leaf surface
    # nodes law - (Ca-Cs) = ra/(ra+rs+rmes)*(Ca-Cc) --> Cs = Ca - ra/(ra+rs+rmes)*(Ca-Cc)
    leaf.Cs = Ca + gleaf*ra/g_m_s_to_mol_m2_s*(Cc-Ca)
    leaf.gleaf = 1 / (ra/g_m_s_to_mol_m2_s + 1.6/leaf.gs + 1/gm)
    leaf.VPD       = max(esat-e_air,1.0); # can be negative at spin up
    leaf.RH        = min(max(e_air/esat,0.001),0.999);    # will need to be corrected alter to define surface RH
    
    #@show leaf.gleaf, Cc

    ΔCc = Ca-leaf.An/leaf.gleaf - Cc
    # Update Cc now
    leaf.Cc = max(1,Ca-leaf.An/leaf.gleaf)
    #leaf.Cc = cinew
    # CiFunc returns the difference between the current Ci and the new Ci
    return ΔCc
end

"""
  setRoughness!(leaf::leaf_params)

Computes leaf roughness lengths, given
- `leaf` One [`leaf_params`](@ref) structure 
"""
function setRoughness!(leaf::leaf_params)
    # adjust roughness coefficients
    leaf.z0m         = 0.1*leaf.height;                     # tree roughness (m)
    leaf.z0h         = 0.1*leaf.height;                     # tree roughness (m) - TODO should be changed later
    leaf.d           = 2/3*leaf.height;                     # tree displacement height (m)
end


"""
  set_leaf_temperature!(mods::PhotoMods, l::leaf_params)

Computes Vcmax, Jmax, Rd, Kc, Ko, esat, desat, at given leaf temperature (if temperature changed), given
- `model` One [`PhotoMods`](@ref) model setup structure
- `leaf` One [`leaf_params`](@ref) structure 
"""
function set_leaf_temperature!(mod::AbstractPhotosynthesis, l::leaf_params)
    #if l.T != l.T_old
      # Adjust Vcmax to leaf T:
      max_carboxylation_rate!(mod.Vmax, l)
      # Adjust Jmax to leaf T:
      max_electron_transport_rate!(mod.Jmax, l)
      # Adjust Michaelis Menten constants:
      michaelis_menten_constants!(mod.MichaelisMenten, l)
      # Respiration rates:
      leaf_respiration!(mod.respiration, l)
      (l.esat, l.desat) = SatVap(l.T);
      l.T_old = l.T
    #end
    # l.kd = max(0.8738,  0.0301*(l.T-273.15)+ 0.0773); # Can implement that later.
end


















#=

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

=#