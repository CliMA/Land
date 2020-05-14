

"Tolerance threshold for Cc iterations"
tol = 0.1
vpd_min = 0.1

" Just a placeholder for now"

Base.@kwdef mutable struct meteo{TT<:Number}
     S_down::TT = -999.;
     L_down::TT = -999.;
     T_air::TT  = -999.;      # T in K
     e_air::TT  = -999.;
     P_air::TT  =  1e5 ;      # surface pressure (Pa)
     Ca::TT     =  400.;
     PAR::TT    = -999.;
     U::TT      = 1e-6;
     zscreen::TT= 10.0; # measurement height - default
     L::TT      = 1e6;  # atmospheric Obukhov length
     # parameter to define stability function for stable case
     stab_type_stable::TT = 1; # 2 Webb correction tends to be unstable at night - suggest not using
     ustar::TT = 1e-6
     g_m_s_to_mol_m2_s::TT = -Inf
     ppm_to_Pa::TT = 0.1
     ra::TT    = 1e6
     APAR::TT = 500.0
    Cs::TT = 0.0
end


Base.@kwdef mutable struct fluxes{TT<:Number}
  Je::TT = 1100.0
  Ac::TT = 0.0
  Aj::TT = 0.0
  Ai::TT = 0.0
  Ap::TT = 0.0
  Ag::TT = 0.0

  Rd::TT = 0.0
  Je_pot::TT = 0.0
  Ja::TT = 0.0
  Je_red::TT = 0.0
  φ::TT = 0.0
  Rn::TT = 0.0
  H::TT = 0.0
  LE::TT = 0.0
  Sap::TT = 0.0
  An_biochemistry::TT = 0.0
  An_diffusion::TT = 0.0

end



# Ball-Berry stomatal conductance model:
function BallBerry!(flux::fluxes, l::leaf_params)
  #  Cs  : CO2 at leaf surface [ppm]
  #  RH  : relative humidity [0-1]
  #  An   : Net assimilation in 'same units of CO2 as Cs' micromoles/m2/s
  #  gs   : moles/m2/s

  l.gs = l.g1_BB * max(flux.An_biochemistry,1e-9) * l.RH/flux.Cs  + l.g0;
  #println("gs=",l.gs,", Cs=",flux.Cs," An_biochemistry=", flux.An_biochemistry," RH=", l.RH )
  if(l.gs<0)
    println("Error - gs=",l.gs,", Cs=",flux.Cs," An_biochemistry=", flux.An_biochemistry," RH=", l.RH )
  end
end # function



# Medlyn stomatal conductance model:
function Medlyn!(flux::fluxes, l::leaf_params)
  #  Cs  : CO2 at leaf surface
  #  VPD  : vapor pressure deficit - Pa
  #  Cs  : CO2 at leaf surface [ppm]
  #  RH  : relative humidity [0-1]
  #  An   : Net assimilation in 'same units of CO2 as Cs' micromoles/m2/s
  #  gs   : moles/m2/s

  l.gs = (1 +l.g1_Medlyn/sqrt(l.VPD)) * max(flux.An_biochemistry,1e-9) /flux.Cs  + l.g0;
end # function


# Gentine stomatal conductance model:
function Gentine!(flux::fluxes, l::leaf_params)
  #  Cs  : CO2 at leaf surface
  #  VPD  : vapor pressure deficit - Pa
  #  Cs  : CO2 at leaf surface [ppm]
  #  RH  : relative humidity [0-1]
  #  An   : Net assimilation in 'same units of CO2 as Cs' micromoles/m2/s
  #  gs   : moles/m2/s

  setLeafkl!(l, l.psi_l) # set hydraulic conductivity of leaf
  l.gs = l.g1_BB*l.kleaf/l.kmax * max(flux.An_biochemistry,1e-9) /flux.Cs  + l.g0;
end # function








"""
    LeafPhotosynthesis!(flux::fluxes, leaf::leaf_params, met::meteo)

Compute net assimilation rate A, fluorescence F using biochemical model

# Arguments
- `flux::fluxes`: fluxes structure.
- `leaf::leaf_params`: leaf_params structure.
- `T::Number`: Leaf Temperature
"""
function LeafPhotosynthesis!(flux::fluxes, leaf::leaf_params, met::meteo)
    # Adjust rates to leaf Temperature (C3 only for now):
    # This is like 1µs already now (used to be 5), need to maybe only run if we change T?
    setLeafT!(leaf)

    # conversion factor for conductance - T and P dependent
    met.g_m_s_to_mol_m2_s = met.P_air/(physcon.Rgas*met.T_air) ; # Körner, C., Scheel, J.A., Bauer, H., 1979. Maximum leaf diffusive conductance in vascular plants. Photosynthetica 13, 45–82.

    # Compute max PSII efficiency here (can later be used with a variable Kn!)
    leaf.Kp = 4.0
    φ_PSII  = leaf.Kp/(leaf.Kp+leaf.Kf+leaf.Kd+leaf.Kn)

    # Save leaf respiration
    flux.Rd = leaf.Rdleaf;

    # Calculate potential electron transport rate (assuming no upper bound, proportional to absorbed light!):
    flux.Je_pot = 0.5 * leaf.maxPSII * flux.APAR;                          # potential electron transport rate (important for later)
    flux.Je_red = 0.5 * φ_PSII * flux.APAR;                                # Includes Kn here
    # Some bound constraint on VPD:
    #flux.ceair = min(max(flux.eair, 0.03*leaf.esat), leaf.esat )

    # Electron transport rate for C3 plants
    # Actual colimited potential Je (curvature and Jmax)
    #@show leaf.θ_j
    #@show flux.Je_red + leaf.Jmax
    #@show flux.Je_red * leaf.Jmax
    flux.Je = minimum(quadratic(leaf.θ_j, -(flux.Je_red + leaf.Jmax), flux.Je_red * leaf.Jmax))    # Bonan eq. 11.21

    if leaf.dynamic_state
      gs_actual = leaf.gs;
      if (leaf.gstyp == 1)
        BallBerry!(flux, leaf)
        #@show leaf.gs
      elseif(leaf.gstyp == 2) # Medlyn default
        Medlyn!(flux, leaf)
      end
      leaf.gs_ss = leaf.gs
      leaf.gs = gs_actual
      #@show leaf.gs
      #@show leaf.gs_ss
      leaf.Cc = CcFuncGs!(gs_actual, flux,leaf, met)
    else
      # computes assimilation
      CcFunc!(flux, leaf, met)
    end

    # Rate of actual CO2 per electron, incl. photorespiration
    # (Cc-Gamma_star)./(Cc+2*Gamma_star)
    leaf.CO2_per_electron = (leaf.Cc-leaf.Γstar)/(leaf.Cc+2leaf.Γstar) * leaf.effcon;

    # Actual effective ETR:
    flux.Ja = max(0,flux.Ag / leaf.CO2_per_electron);
    flux.Ja = min(flux.Ja,flux.Je_pot )

    # Effective photochemical yield:
    flux.φ = leaf.maxPSII*flux.Ja/flux.Je_pot;
    #println(flux.Ja, " ", flux.Je_pot)
    flux.φ = min(1/leaf.maxPSII,flux.φ)
    x   = max(0,  1-flux.φ/leaf.maxPSII);       # degree of light saturation: 'x' (van der Tol e.Ap. 2014)
    Fluorescencemodel!(flux.φ,x,leaf)
    #@show leaf.gs

end # LeafPhotosynthesis (similar to biochem in SCOPE)





function setRoughness!(leaf::leaf_params)
    # adjust roughness coefficients
    leaf.z0m         = 0.1*leaf.height;                     # tree roughness (m)
    leaf.z0h         = 0.1*leaf.height;                     # tree roughness (m) - TODO should be changed later
    leaf.d           = 2/3*leaf.height;                     # tree displacement height (m)
end



"""
    CcFunc!(flux::fluxes, leaf::leaf_params, met::meteo)

Compute Assimilation using Cc as input

# Arguments
- `Cc::Number`: Cc.
- `flux::fluxes`: fluxes structure.
- `leaf::leaf_params`: leaf_params structure.
"""
function CcFunc!(flux::fluxes, leaf::leaf_params, met::meteo)
    @unpack  Γstar, Cc, Kc, Ko, o₂, Vcmax = leaf
    #println("Cc_in=",leaf.Cc)
    if leaf.C3
        # C3: Rubisco-limited photosynthesis; still need to check CO2 mixing ratios vs partial pressures.
        # still need to include ppm2bar (but can be done on leaf structure!)
        flux.Ac = Vcmax * max(Cc-Γstar, 0) / (Cc + Kc*(1 +o₂/Ko)) # Bonan eq. 11.28
        # C3: RuBP-limited photosynthesis (this is the NADPH requirement stochiometry)
        flux.Aj = flux.Je * max(Cc-Γstar, 0) / (4Cc + 8Γstar)     # Bonan eq. 11.29

        # for C3, set ap to Inf
        flux.Ap = Inf
    else #C4 Photosynthesis, still to be implemented
        flux.Ac = flux.Aj = flux.Ap = 0
    end
    # Net photosynthesis as the minimum or co-limited rate
    if leaf.use_colim
        flux.Ai = minimum(quadratic(leaf.C3 ? 0.98 : 0.80, -(flux.Ac + flux.Aj), flux.Ac * flux.Aj))
        if leaf.C3
            flux.Ag = flux.Ai
        else # C4 colimitation with ap
            flux.Ag = minimum(quadratic(0.95, -(flux.Ai + flux.Ap), flux.Ai*flux.Ap))                 # Bonan Eq 11.33
        end
    else
        flux.Ag = min(flux.Ac,flux.Aj,flux.Ap) # gross assimilation
    end
    # Prevent photosynthesis from ever being negative
    flux.Ag = max(0,flux.Ag)
    flux.Ai = max(0,flux.Ai)
    flux.Aj = max(0,flux.Aj)
    flux.Ap = max(0,flux.Ap)

    # Net photosynthesis due to biochemistry
    flux.An_biochemistry = flux.Ag - leaf.Rdleaf # net assimilation

    # adjust aerodynamic resistance based on leaf boundary layer and Monin Obukhov
    setra!(leaf, flux, met)

    # CO2 at leaf surface
    # nodes law - (Ca-Cs) = ra/(ra+rs+rmes)*(Ca-Cc) --> Cs = Ca - ra/(ra+rs+rmes)*(Ca-Cc)
    leaf.gleaf = 1 / (flux.ra/flux.g_m_s_to_mol_m2_s + 1.6/leaf.gs + 1.0/leaf.gm)
    flux.Cs = met.Ca + leaf.gleaf*flux.ra/flux.g_m_s_to_mol_m2_s*(leaf.Cc-met.Ca)

    # println("Cs=",flux.Cs,", ra=",flux.ra, ", ra/rleaf=",leaf.gleaf*flux.ra/flux.g_m_s_to_mol_m2_s, ", Cc=", leaf.Cc, ", Ca=", met.Ca, " L=", met.L, " u*=", flux.ustar, " H=",flux.H)

    # compute stomatal conductance gs
    leaf.VPD       = max(leaf.esat-met.e_air,1.0); # can be negative at spin up
    leaf.RH        = min(max(met.e_air/leaf.esat,0.001),0.999);    # will need to be corrected alter to define surface RH

    if (leaf.gstyp == 1)
        BallBerry!(flux, leaf)
    elseif(leaf.gstyp == 2) # Medlyn default
        Medlyn!(flux, leaf)
    elseif(leaf.gstyp == 3) # Medlyn default
        Gentine!(flux, leaf)
    end

    # Diffusion (supply-based) photosynthetic rate - Calculate Cc from the diffusion rate
    # total conductance - mol/m2/s

    flux.An_diffusion = leaf.gleaf*(met.Ca - leaf.Cc)

end

"""
    CcFuncGs!(gs::Number, flux::fluxes, leaf::leaf_params)

Compute Assimilation using fixed <stomatal></stomatal> conductance gs.
# Arguments
- `gs::Number`: Stomatal conductance.
- `flux::fluxes`: fluxes structure.
- `leaf::leaf_params`: leaf_params structure.
"""
function CcFuncGs!(gs::Number, flux::fluxes, leaf::leaf_params, met::meteo)
    FT = eltype(gs)
    @unpack  Γstar, Cc, Kc, Ko, o₂, Vcmax, Rdleaf = leaf
    @unpack Ca = met
    # Compute overall conductance (Boundary layer, stomata and mesophyll)
    gleaf = 1 / (flux.ra/flux.g_m_s_to_mol_m2_s + FT(1.6)/gs + FT(1.0)/leaf.gm)
    if gleaf<eps(FT) gleaf=eps(FT) end
    leaf.gleaf = gleaf
    #flux.ac = leaf.vcmax * max(Ci-leaf.Γstar, 0.0) / (Ci + leaf.kc*(1.0+leaf.o₂/leaf.ko)) # Bonan eq. 11.28
    # C3: RuBP-limited photosynthesis (this is the NADPH requirement stochiometry)
    #flux.aj = flux.je * max(Ci-leaf.Γstar, 0.0) / (4.0*Ci + 8.0*leaf.Γstar)               # Bonan eq. 11.29
    if leaf.C3
        # C3: Rubisco-limited photosynthesis; still need to check CO2 mixing ratios vs partial pressures.
        # still need to include ppm2bar (but can be done on leaf structure!)
        flux.Ac = Vcmax * max(Cc-Γstar, 0) / (Cc + Kc*(1 +o₂/Ko)) # Bonan eq. 11.28
        # C3: RuBP-limited photosynthesis (this is the NADPH requirement stochiometry)
        flux.Aj = flux.Je * max(Cc-Γstar, 0) / (4Cc + 8Γstar)     # Bonan eq. 11.29
        # Photosynthesis limited by diffusion:
        #flux.An_diffusion = 10.0;#gleaf*(Ca - Cc)
        # for C3, set ap to Inf
        flux.Ap = Inf
    else #C4 Photosynthesis, still to be implemented
        flux.Ac = flux.Aj = flux.Ap = 0
    end

    flux.Ag = min(flux.Ac,flux.Aj,flux.Ap)
    
    flux.Ag = max(0,flux.Ag)
    flux.Ai = max(0,flux.Ai)
    flux.Aj = max(0,flux.Aj)
    flux.Ap = max(0,flux.Ap)
    
    # Compute net Photosynthesis
    flux.An_biochemistry = flux.Ag - Rdleaf
    # adjust aerodynamic resistance based on leaf boundary layer and Monin Obukhov
    setra!(leaf, flux, met)
    #@show flux.LE
    # CO2 at leaf surface
    # nodes law - (Ca-Cs) = ra/(ra+rs+rmes)*(Ca-Cc) --> Cs = Ca - ra/(ra+rs+rmes)*(Ca-Cc)
    leaf.gleaf = 1 / (flux.ra/flux.g_m_s_to_mol_m2_s + FT(1.6)/leaf.gs + FT(1.0)/leaf.gm)
    flux.Cs = met.Ca + leaf.gleaf*flux.ra/flux.g_m_s_to_mol_m2_s*(leaf.Cc-met.Ca)

    # println("Cs=",flux.Cs,", ra=",flux.ra, ", ra/rleaf=",leaf.gleaf*flux.ra/flux.g_m_s_to_mol_m2_s, ", Cc=", leaf.Cc, ", Ca=", met.Ca, " L=", met.L, " u*=", flux.ustar, " H=",flux.H)

    # compute stomatal conductance gs
    leaf.VPD       = max(leaf.esat-met.e_air,1.0); # can be negative at spin up
    leaf.RH        = min(max(met.e_air/leaf.esat,0.001),0.999);    # will need to be corrected alter to define surface RH

    # Compute CO2 at leaf surface
    flux.Cs = Ca - flux.An_biochemistry / (flux.g_m_s_to_mol_m2_s/flux.ra)

    # Compute Cc (included Mesophyll as well in principle)
    Cc_val = Ca - flux.An_biochemistry / gleaf
    #leaf.CO2_per_electron = (ci_val-leaf.Γstar)./(ci_val+2.0*leaf.Γstar) .* leaf.effcon;
end # Function CiFuncGs!


function ψ_m(ζ,stab_type_stable) # momentum correction function
  FT = eltype(ζ)

    if(ζ>=0) # stable conditions
        if(stab_type_stable==1) # Grachev default value - Grachev et al SHEBA 2007
          x  = (1+ζ)^FT(1/3);
          am = FT(5.0);
          bm = am/FT(6.5);
          ah = bh = FT(5.0);
          ch = FT(3.0);
          Bm = ((FT(1.0)-bm)/bm)^FT(1/3);
          Bh = sqrt(5);
          ψ = -FT(3)*am/bm*(x-1) + am*Bm/(2bm)*( FT(2)*log((x+Bm)/(1+Bm))
                                                   - log( (x*x-x*Bm+Bm*Bm)/(1-Bm+Bm*Bm) )
                                                   +2*sqrt(FT(3.0))*( atan((2x-Bm)/(sqrt(FT(3.0))*Bm)) -  atan((2-Bm)/(sqrt(FT(3.0))*Bm)))  );
        elseif(stab_type_stable==2) # Webb correction as in NoahMP, quite numerically unstable at night - lots of oscillations
          alpha = FT(5.0);
          """ should be the correct one but not continuously differntialble
          if(ζ>=0.0 && ζ<=1)
            ψ = -alpha * ζ;
          else
            ψ = -alpha; # should be z less anyways but here just made to create a plateau in the value
          end
          """
          # Pierre's new one:
          ψ = -alpha*tanh(ζ);
        elseif(stab_type_stable==3) # Beljaars-Holtslag, 1991, Journal of Applied Meteorology
          # Constants
          a = FT(1.0);
          b = FT(0.667);
          c = FT(5);
          d = FT(0.35);
          # Stability correction for momentum
          ψ = -(a * ζ
                   + b * (ζ - c / d) * exp(-d * ζ)
                   + b * c / d);
        elseif(stab_type_stable==4) #Cheng and Brutsaert (2005) as described by Jimenez et al. (2012)
          a = FT(6.1);
          b = FT(2.5);
          c = FT(5.3);
          d = FT(1.1);
          ψ = -a * log(ζ + (1 + ζ^b) ^ (1 / b));
        end

    else #  (ζ<0) unstable
        x   =   (1-16ζ)^FT(0.25);
        ψ   =   FT(2.0)*log((1+x)/2) + log((1+x*x)/2) - 2*atan(x) + FT(π)/2;
    end
end




function ψ_h(ζ,stab_type_stable) # momentum correction function

if(ζ>=0.0) # stable conditions
      if(stab_type_stable==1) # Grachev default value - Grachev et al SHEBA 2007
        # issue when x ~ 0
        x  = (1.0+ζ)^(1.0/3.0);
        am = 5.0;
        bm = am/6.5;
        ah = bh = 5.0;
        ch = 3.0;
        Bm = ((1.0-bm)/bm)^(1.0/3.0);
        Bh = sqrt(5);
        ψ = -bh/2*log(1+ch*ζ+ζ*ζ) + (-ah/Bh + bh*ch/(2.0*Bh))*(log((2.0*ζ+ch-Bh)/(2.0*ζ+ch+Bh)) - log((ch-Bh)/(ch+Bh)) );
      elseif(stab_type_stable==2) # Webb correction as in NoahMP
        alpha = 5.0;
        #if(ζ>=0.0 && ζ<=1)
        #  ψ = -alpha * ζ;
        #else
        #  ψ = -alpha;
        #end
        # Pierre - I didn't like this so I changed it to be C1
        ψ = -alpha*tanh(ζ);
      elseif(stab_type_stable==3) # Beljaars-Holtslag, 1991, Journal of Applied Meteorology
        # Constants
        a = 1.0;
        b = 0.667;
        c = 5;
        d = 0.35;
        # Stability correction for momentum
        ψ = (- (1. + (2. / 3.) * a * ζ) ^ (3. / 2.)
                  - b * (ζ - (c / d)) * exp(-d * ζ)
                  - (b * c) / d + 1);
      elseif(stab_type_stable==4) #Cheng and Brutsaert (2005) as described by Jimenez et al. (2012)
        a = 6.1;
        b = 2.5;
        c = 5.3;
        d = 1.1;
        ψ =  -c * log(ζ + (1 + ζ^d) ^ (1. / d));
      end

    else # (ζ<0) unstable
        x   = (1.0-16.0*ζ)^0.25;
        ψ = 2.0*log((1.0+x*x)/2.0);
    end

    return ψ
end

function setra!(l::leaf_params, flux::fluxes, met::meteo) # set aerodynamic resistance
    # based on Monin-Obukhov Similiarity theory -> to be changed for LES
    # compute Obukhov length
    # iterate a few times

    # TODO ideally should not have any information about the leaves as it is the turbuence above canopy in log profile
    # first update roughness (if phenology is changing)
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
    DeltaT   =   l.T - met.T_air;
    esat,desat_dT = SatVap(l.T);
    VPD      =   max(esat-met.e_air,1.0);
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
      ram_full =   ra_leaf + ra_m;
      raw_full =   ra_leaf + ra_w;
      H        =   ρd*physcon.Cpd*DeltaT/raw_full;
      rs_s_m   =   flux.g_m_s_to_mol_m2_s/l.gs;
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
    flux.H  = H
    flux.LE = LE
    flux.ustar = ustar
    flux.ra = raw_full
end


