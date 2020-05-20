
const FT = Float32
abstract type AbstractPhotosynthesis end

"""
PhotoMods
describes all necessary modules
"""
Base.@kwdef struct PhotoMods{FM,PM,RM,SM,JM,VM,MM,BL,CL} <: AbstractPhotosynthesis
    fluorescence::FM       = FlexasTolBerryFluorescence{FT}()
    photosynthesis::PM     = C3FvCBPhoto()
    respiration::RM        = RespirationCLM{FT}()
    stomatal::SM           = BallBerryStomata{FT}()
    Jmax::JM               = JmaxCLM{FT}()
    Vmax::VM               = VcmaxCLM{FT}()
    MichaelisMenten::MM    = MM_CLM{FT}()
    BoundaryLayer::BL      = GentineLeafBoundary()
    colimitation::CL     = CurvedColimit{FT}()
end

"""
    LeafPhotosynthesis!(mod::AbstractPhotosynthesis, leaf::leaf_params, met::meteo)

Compute net assimilation rate A, fluorescence F using biochemical model, given
- `model` One [`PhotoMods`](@ref) model setup structure
- `leaf` One [`leaf_params`](@ref) structure 
- `met` One [`lmeteo`](@ref) structure 

"""
function LeafPhotosynthesis!(mo::AbstractPhotosynthesis, leaf::leaf_params, met::meteo)
  # Set leaf temperature:
  set_leaf_temperature!(mo, leaf)
  met.g_m_s_to_mol_m2_s = (met.P_air-met.e_air)/(physcon.Rgas*met.T_air) 
  met.ppm_to_Pa = (met.P_air-met.e_air)*1e-6

  
  # Compute Cc and Photosynthesis:
  if leaf.dynamic_state
    CcFunc!(leaf.Cc; mods=mo, leaf=leaf, met=met)
  else
    @inline f(x) = CcFunc!(x; mods=mo, leaf=leaf, met=met)
    sol = find_zero(f, SecantMethod{FT}(0.0, 2met.Ca), CompactSolution(),SolutionTolerance{FT}(1e-3))
    leaf.Cc = sol.root
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
    leaf.Cc = Cc
    # Compute Rubisco Limited Photosynthesis
    rubisco_limited_rate!(mods.photosynthesis, leaf, met)
    # Light limited rate
    light_limited_rate!(mods.photosynthesis, leaf, met, leaf.APAR)
    # Product limited rate
    product_limited_rate!(mods.photosynthesis, leaf, met)
    
    # adjust aerodynamic resistance based on leaf boundary layer and Monin Obukhov
    boundary_layer_resistance!(mods.BoundaryLayer, leaf,  met)

    @unpack Ac, Aj,Ap, Cs, VPD, RH, gm, gs, Ag, Rd, An, gleaf, esat, Cc = leaf
    @unpack Ca, ra, g_m_s_to_mol_m2_s, e_air = met

    # add colimitation later:
    #@show Ac,Aj,Ap, Cc
    leaf.Ag = photosynthesis_colimit(mods.colimitation, Ac, Aj, Ap) 
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




