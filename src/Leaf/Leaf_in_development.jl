###############################################################################
#
# Update leaf photosynthesis
# Move to SPAC because of its dependence on Plant sub-module
#
###############################################################################
# TODO need to abstract this one for Gentine model
function leaf_photosynthesis!(paraset::C3ParaSet, leaf::LeafParams, met::MeteoParams, APAR::FT, model::ESMBallBerry) where{FT}
    if !leaf.dynamic_state
        # need to call Ball-Berry model type in Plant
        # TODO redo this function some time using CO2 directly in the get_empirical_gsw function
        leaf.gs, leaf.p_i = get_empirical_gsw_pi(leaf.θ_j, FT(0.0), FT(0.3), FT(0.8), leaf.Jmax25, leaf.Vpmax25, met.p_a, met.p_atm, met.e_air,leaf.p_O₂, APAR, FT(0.4), leaf.Rd25, leaf.T, leaf.Vcmax25, paraset, model)
        # TODO negelecting ra for now, fix it later!
        leaf.gleaf = 1 / (1.6/leaf.gs + 1/leaf.gm)
    end
    update_leaf_TD!(paraset, leaf)
    electron_transport_rate!(paraset, leaf, APAR)
    rubisco_limited_rate!(paraset, leaf)
    light_limited_rate!(paraset, leaf, APAR)
    product_limited_rate!(paraset, leaf)
    leaf.Ag = min(leaf.Aj, leaf.Ac, leaf.Ap)
    leaf.An = leaf.Ag - leaf.Rd
end

function leaf_photosynthesis!(paraset::C4ParaSet, leaf::LeafParams, met::MeteoParams, APAR::FT, model::ESMBallBerry) where{FT}
    if !leaf.dynamic_state
        # need to call Ball-Berry model type in Plant
        # TODO redo this function some time
        leaf.gs, leaf.p_i = get_empirical_gsw_pi(leaf.θ_j, FT(0.0), FT(0.3), FT(0.8), leaf.Jmax25, leaf.Vpmax25, met.p_a, met.p_atm, met.e_air,leaf.p_O₂, APAR, FT(0.4), leaf.Rd25, leaf.T, leaf.Vcmax25, paraset, model)
        # TODO negelecting ra for now, fix it later!
        leaf.gleaf = 1 / (1.6/leaf.gs + 1/leaf.gm)
    end
    update_leaf_TD!(paraset, leaf)
    rubisco_limited_rate!(paraset, leaf)
    light_limited_rate!(paraset, leaf, APAR)
    product_limited_rate!(paraset, leaf)
    leaf.Ag = min(leaf.Aj, leaf.Ac, leaf.Ap)
    leaf.An = leaf.Ag - leaf.Rd
end





###############################################################################
#
# Update the boundary layer conductance of leaf
#
###############################################################################
"""
    boundary_layer_resistance!(bl::LeafBLParaSetFixed, l::LeafParams,  met::meteo)

Apply fixed leaf boundary layer resistance, given
- `bl` A `LeafBLParaSetFixed` or `LeafBLParaSetGentine` type leaf boundary layer
- `leaf` A `LeafParams` struct
- `met` A `meteo` type parameter set for temporary variables
"""
function boundary_layer_resistance!(bl::LeafBLParaSetFixed, leaf::LeafParams,  met::MeteoParams)
    # compute the evolving buoyancy flux
    Tv       = met.T_air * (1.0 + 0.61*met.e_air)
    ΔT       = leaf.T - met.T_air
    VPD      = max(leaf.esat-met.e_air, 1.0)
    ρd       = met.P_air / (R_D * met.T_air)

    # TODO use ThermoDynamics to compute latent heat flux
    lv       = latent_heat_vapor(leaf.T)
    _ra      = bl.ra
    leaf.H      = ρd * CP_D * ΔT / _ra;

    # TODO Fix this?
    rs_s_m   = met.g_m_s_to_mol_m2_s/leaf.gs;
    leaf.LE     = WATER_AIR_MRATIO / met.P_air * ρd * lv * VPD / (rs_s_m + _ra)

    # TODO  _ra is constant, no update is necessary
    met.ra = _ra
end

function boundary_layer_resistance!(bl::LeafBLParaSetGentine, leaf::LeafParams,  met::MeteoParams)
    # based on Monin-Obukhov Similiarity theory -> to be changed for LES
    # compute Obukhov length
    # iterate a few times

    # TODO ideally should not have any information about the leaves as it is the turbuence above canopy in log profile
    # first update roughness (if phenology is changing)
    # need to make this abstract as well.
    setRoughness!(leaf)

    if(leaf.height>met.zscreen)
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
    leaf.Cd     = Inf;


    # need to compute the evolving buoyancy flux
    Tv       = met.T_air*(1.0+0.61*met.e_air);
    ra_leaf  = 1/leaf.Cd / sqrt(met.U/leaf.dleaf);
    #println("ra_leaf=",ra_leaf)
    ΔT   =   leaf.T - met.T_air;
    VPD      =   max(leaf.esat-met.e_air,1.0)
    ρd       =   met.P_air/(R_D * met.T_air)

    # TODO use ThermoDynamics to compute latent heat flux
    lv       =   latent_heat_vapor(leaf.T);
    L        =   met.L; # initial Obukhov length
    counter = 1;
    #@show abs(1.0-Lold/L)
    while (counter<20 && abs(1.0-Lold/L)>1e-4) # 1% error
      #println("L=",L," ,Lold=",Lold)
      ra_m     =   max(1.0/(VON_KARMAN_CONST^2*met.U) * ( log((met.zscreen - leaf.d)/leaf.z0m) - ψ_m((met.zscreen - leaf.d)/L,met.stab_type_stable) + ψ_m(leaf.z0m/L,met.stab_type_stable) ) * ( log((met.zscreen - leaf.d)/leaf.z0m) - ψ_m((met.zscreen - leaf.d)/L,met.stab_type_stable) + ψ_h(leaf.z0h/L,met.stab_type_stable) ), rmin) ;# momentum aerodynamic resistance
      ra_w     =   max(1.0/(VON_KARMAN_CONST^2*met.U) * ( log((met.zscreen - leaf.d)/leaf.z0m) - ψ_m((met.zscreen - leaf.d)/L,met.stab_type_stable) + ψ_m(leaf.z0m/L,met.stab_type_stable) ) * ( log((met.zscreen - leaf.d)/leaf.z0h) - ψ_h((met.zscreen - leaf.d)/L,met.stab_type_stable) + ψ_h(leaf.z0h/L,met.stab_type_stable) ), rmin) ;# water aerodynamic resistance
      #println("ra_m=",ra_m)
      #@show ra_leaf
      #@show ra_w
      ram_full =   ra_leaf + ra_m;
      raw_full =   ra_leaf + ra_w;
      H        =   ρd * CP_D * ΔT / raw_full;
      rs_s_m   =   met.g_m_s_to_mol_m2_s/leaf.gs;
      LE       =   WATER_AIR_MRATIO / met.P_air * ρd * lv * VPD / (rs_s_m + raw_full)
      #@show LE
      ustar    =   sqrt(met.U/ram_full);
      Hv_s     =   H + 0.61 * CP_D / lv * met.T_air * LE
      Lold     =   L;
      L        =   - ustar^3*Tv/(GRAVITY * VON_KARMAN_CONST * Hv_s); # update Obukhov length
      #println("L=",L, " ra=",ra_w," (s/m), H=", H, " (W/m2), counter=", counter)
      counter = counter+1
    end

    # TODO Is the met more a "global" variable set for the tree? Or a "local" variable for each layer/leaf?
    # If for the whole plant or layer, then met should be updated only once?
    # save these values in leaf and flux structures
    met.L   = L
    leaf.H  = H
    leaf.LE = LE
    met.ustar = ustar
    met.ra = raw_full
end







"""
ψ_h(ζ,stab_type_stable)

momentum correction function
"""
function ψ_h(ζ,stab_type_stable)
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




"""
ψ_m(ζ,stab_type_stable)

Momentum correction function
"""
function ψ_m(ζ,stab_type_stable)
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

return ψ
end




"""
setRoughness!(leaf::LeafPhotosynthesisParams)

Computes leaf roughness lengths, given
- `leaf` One [`LeafPhotosynthesisParams`](@ref) structure 
"""
function setRoughness!(leaf::LeafParams)
# adjust roughness coefficients
leaf.z0m         = 0.1*leaf.height;                     # tree roughness (m)
leaf.z0h         = 0.1*leaf.height;                     # tree roughness (m) - TODO should be changed later
leaf.d           = 2/3*leaf.height;                     # tree displacement height (m)
end
