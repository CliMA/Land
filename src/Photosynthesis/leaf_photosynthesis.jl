
const FT = Float32
abstract type AbstractPhotosynthesis end
abstract type AbstractLeafBoundaryLayer end
struct GentineLeafBoundary <: AbstractLeafBoundaryLayer end

Base.@kwdef struct PhotoMods{FM,PM,RM,SM,JM,VM,MM,BL} <: AbstractPhotosynthesis
    fluorescence::FM       = FlexasTolBerryFluorescence{FT}()
    photosynthesis::PM     = C3FvCBPhoto()
    respiration::RM        = RespirationCLM{FT}()
    stomatal::SM           = BallBerryStomata{FT}()
    Jmax::JM               = JmaxCLM{FT}()
    Vmax::VM               = VcmaxCLM{FT}()
    MichaelisMenten::MM    = MM_CLM{FT}()
    BoundaryLayer::BL      = GentineLeafBoundary()
end

"""
    LeafPhotosynthesis!(mod::AbstractPhotosynthesis, leaf::leaf_params, met::meteo)

Compute net assimilation rate A, fluorescence F using biochemical model

"""
function LeafPhotosynthesis!(mo::AbstractPhotosynthesis, leaf::leaf_params, met::meteo)
  # Set leaf temperature:
  set_leaf_temperature!(mo, leaf)
  met.g_m_s_to_mol_m2_s = met.P_air/(physcon.Rgas*met.T_air) 
  met.ppm_to_Pa = met.P_air*1e-6

  # Compute Cc and Photosynthesis:
  CcFunc!(mo, leaf, met)
  leaf.Cc = met.Ca-leaf.An/leaf.gleaf
  
  # Compute Fluorescence
  leaf_fluorescence!(mo.fluorescence,  leaf)

end # LeafPhotosynthesis (similar to biochem in SCOPE)


"""
    CcFunc!(mod::AbstractPhotosynthesis,  leaf::leaf_params, met::meteo)

Compute Assimilation using Cc as input
"""
function CcFunc!(mods::AbstractPhotosynthesis,  leaf::leaf_params, met::meteo)
    @unpack Ac,Aj,Ap = leaf
    # Compute Rubisco Limited Photosynthesis
    rubisco_limited_rate!(mods.photosynthesis, leaf, met)
    # Light limited rate
    light_limited_rate!(mods.photosynthesis, leaf, met, leaf.APAR)
    # Product limited rate
    product_limited_rate!(mods.photosynthesis, leaf)
    
    # add colimitation later:
    leaf.Ag = min(Ac,Aj,Ap) # gross assimilation
    #@show leaf.Ag
    # Net photosynthesis due to biochemistry
    leaf.An = leaf.Ag - leaf.Rd # net assimilation

    # adjust aerodynamic resistance based on leaf boundary layer and Monin Obukhov
    boundary_layer_resistance!(mods.BoundaryLayer, leaf,  met)

    # CO2 at leaf surface
    # nodes law - (Ca-Cs) = ra/(ra+rs+rmes)*(Ca-Cc) --> Cs = Ca - ra/(ra+rs+rmes)*(Ca-Cc)
    leaf.gleaf = 1 / (met.ra/met.g_m_s_to_mol_m2_s + 1.6/leaf.gs + 1.0/leaf.gm)
    leaf.Cs = met.Ca + leaf.gleaf*met.ra/met.g_m_s_to_mol_m2_s*(leaf.Cc-met.Ca)

    # println("Cs=",flux.Cs,", ra=",flux.ra, ", ra/rleaf=",leaf.gleaf*flux.ra/flux.g_m_s_to_mol_m2_s, ", Cc=", leaf.Cc, ", Ca=", met.Ca, " L=", met.L, " u*=", flux.ustar, " H=",flux.H)

    # compute stomatal conductance gs
    leaf.VPD       = max(leaf.esat-met.e_air,1.0); # can be negative at spin up
    leaf.RH        = min(max(met.e_air/leaf.esat,0.001),0.999);    # will need to be corrected alter to define surface RH

    # Compute stomatal conductance:
    stomatal_conductance!(mods.stomatal, leaf)

end


function setRoughness!(leaf::leaf_params)
    # adjust roughness coefficients
    leaf.z0m         = 0.1*leaf.height;                     # tree roughness (m)
    leaf.z0h         = 0.1*leaf.height;                     # tree roughness (m) - TODO should be changed later
    leaf.d           = 2/3*leaf.height;                     # tree displacement height (m)
end








# Set Leaf rates with Vcmax, Jmax and rd at 25C as well as actual T here:
# For some reason, this is slow and allocates a lot, can be improved!!
"Set Leaf rates with Vcmax, Jmax and rd at 25C as well as actual T here"
function set_leaf_temperature!(mod::AbstractPhotosynthesis, l::leaf_params)
    if l.T != l.T_old
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
    end
    # l.kd = max(0.8738,  0.0301*(l.T-273.15)+ 0.0773); # Can implement that later.
end

function boundary_layer_resistance!(mod::GentineLeafBoundary, l::leaf_params,  met::meteo) # set aerodynamic resistance
    # based on Monin-Obukhov Similiarity theory -> to be changed for LES
    # compute Obukhov length
    # iterate a few times

    # TODO ideally should not have any information about the leaves as it is the turbuence above canopy in log profile
    # first update roughness (if phenology is changing)
    # need to make this abstract as well.
    setRoughness!(l)

    if(l.height>met.zscreen)
        println("Vegetation height is higher than screen level height")
        #process.exit(20)
    end

    rmin  = 10.0;
    Lold  = -1e6;
    raw_full = -999.0;
    H     = -999.0;
    LE    = -999.0;
    ustar = -999.0;

    # TODO change this later - 0 leaf boundary layer resistance - I prefer this as it is somewhat included in z0
    l.Cd     = Inf;


    # need to compute the evolving buoyancy flux
    Tv       = met.T_air*(1.0+0.61*met.e_air);
    ra_leaf  = 1/l.Cd / sqrt(met.U/l.dleaf);
    #println("ra_leaf=",ra_leaf)
    ΔT   =   l.T - met.T_air;
    VPD      =   max(l.esat-met.e_air,1.0);
    #@show VPD
    ρd       =   met.P_air/(physcon.Rd*met.T_air);      # dry air density (kg/m3)
    lv       =   Lv(l.T);
    L        =   met.L; # initial Obukhov length
    counter = 1;
    #@show abs(1.0-Lold/L)
    while (counter<20 && abs(1.0-Lold/L)>1e-4) # 1% error
      #println("L=",L," ,Lold=",Lold)
      ra_m     =   max(1.0/(physcon.K^2*met.U) * ( log((met.zscreen - l.d)/l.z0m) - ψ_m((met.zscreen - l.d)/L,met.stab_type_stable) + ψ_m(l.z0m/L,met.stab_type_stable) ) * ( log((met.zscreen - l.d)/l.z0m) - ψ_m((met.zscreen - l.d)/L,met.stab_type_stable) + ψ_h(l.z0h/L,met.stab_type_stable) ), rmin) ;# momentum aerodynamic resistance
      ra_w     =   max(1.0/(physcon.K^2*met.U) * ( log((met.zscreen - l.d)/l.z0m) - ψ_m((met.zscreen - l.d)/L,met.stab_type_stable) + ψ_m(l.z0m/L,met.stab_type_stable) ) * ( log((met.zscreen - l.d)/l.z0h) - ψ_h((met.zscreen - l.d)/L,met.stab_type_stable) + ψ_h(l.z0h/L,met.stab_type_stable) ), rmin) ;# water aerodynamic resistance
      #println("ra_m=",ra_m)
      #@show ra_leaf
      #@show ra_w
      ram_full =   ra_leaf + ra_m;
      raw_full =   ra_leaf + ra_w;
      H        =   ρd*physcon.Cpd*ΔT/raw_full;
      rs_s_m   =   met.g_m_s_to_mol_m2_s/l.gs;
      LE       =   physcon.ε/met.P_air*ρd*lv*VPD/(rs_s_m+raw_full);
      #@show LE
      ustar    =   sqrt(met.U/ram_full);
      Hv_s     =   H + 0.61 * physcon.Cpd/lv * met.T_air * LE;
      Lold     =   L;
      L        =   - ustar^3*Tv/(physcon.grav*physcon.K*Hv_s); # update Obukhov length
      #println("L=",L, " ra=",ra_w," (s/m), H=", H, " (W/m2), counter=", counter)
      counter = counter+1
    end

    # save these values in leaf and flux structures
    met.L   = L
    l.H  = H
    l.LE = LE
    met.ustar = ustar
    met.ra = raw_full
end


