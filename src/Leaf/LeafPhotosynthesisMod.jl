module LeafPhotosynthesisMod

using Parameters

using ..PhysCon
using ..WaterVaporMod

include("Leaf.jl")
include("surface_energy_balance.jl")

export fluxes
"Tolerance thhreshold for Ci iterations"
tol = 0.1
vpd_min = 0.1

" Just a placeholder for now"
@with_kw mutable struct fluxes{TT<:Number}
         APAR::TT = 500.0
         gbc::TT = 100.0
         gbv::TT = 100.0
         ceair::TT = 1400.0
         eair::TT = 1400.0
         Je::TT = 1100.0
         Ac::TT = 0.0
         Aj::TT = 0.0
         Ai::TT = 0.0
         Ap::TT = 0.0
         Ag::TT = 0.0
         Cs::TT = 0.0
         #ci::TT = 0.0
         rd::TT = 0.0
         Je_pot::TT = 0.0
         Ja::TT = 0.0
         Je_red::TT = 0.0
         φ::TT = 0.0
         Rn::TT = 0.0
         H::TT = 0.0
         LE::TT = 0.0
         Sap::TT = 0.0
         An::TT = 0.0
         ustar::TT = 1e-6
         g_m_s_to_micromol_m2_s::TT = 1.0/38.0e6
end



# Ball-Berry stomatal conductance model:
function BallBerry!(flux::fluxes,  l::leaf_params)
  #  Cs  : CO2 at leaf surface
  #  RH  : relative humidity
  #  A   : Net assimilation in 'same units of CO2 as Cs'/m2/s
  #  minCi : minimum Ci as a fraction of Cs (in case RH is very low?)
  #  Ci_input : will calculate gs if A is specified.

  l.gs = l.g1 * max(flux.An,0.0) * l.RH/flux.Cs  + l.g0;
end # function



# Medlyn stomatal conductance model:
function Medlyn!(flux::fluxes, l::leaf_params)
  #  Cs  : CO2 at leaf surface
  #  VPD  : vapor pressure deficit - Medlyn model
  #  A   : Net assimilation in 'same units of CO2 as Cs'/m2/s
  #  minCi : minimum Ci as a fraction of Cs (in case RH is very low?)
  #  Ci_input : will calculate gs if A is specified.

  #print(A,Cs,l.g0,l.g1,l.gs)
  l.gs = (1+l.g1/sqrt(l.VPD)) * max(flux.An,0.0) /flux.Cs  + l.g0;
  #print(l.gs)
end # function







"""
    LeafPhotosynthesis(flux::fluxes, leaf::leaf_params, met::meteo)

Compute net assimilation rate A, fluorescence F using biochemical model

# Arguments
- `flux::fluxes`: fluxes structure.
- `leaf::leaf_params`: leaf_params structure.
- `T::Number`: Leaf Temperature
"""
function LeafPhotosynthesis(flux::fluxes, leaf::leaf_params, met::meteo)
    # Adjust rates to leaf Temperature (C3 only for now):
    setLeafT!(leaf)

    flux.g_m_s_to_micromol_m2_s = 1.0e6/(44.6*273.15/met.T_air*met.P_air/1.013e5); # Körner, C., Scheel, J.A., Bauer, H., 1979. Maximum leaf diffusive conductance in vascular plants. Photosynthetica 13, 45–82.

    # Compute max PSII efficiency here (can later be used with a variable Kn!)
    leaf.Kp = 4.0
    φ_PSII  = leaf.Kp/(leaf.Kp+leaf.Kf+leaf.Kd+leaf.Kn)

    # Save leaf respiration
    flux.rd = leaf.Rdleaf;
    # Calculate potential electron transport rate (assuming no upper bound, proportional to absorbed light!):
    flux.Je_pot = 0.5 * leaf.maxPSII * flux.APAR;                          # potential electron transport rate (important for later)
    flux.Je_red = 0.5 * φ_PSII * flux.APAR;                                # Includes Kn here
    # Some bound constraint on VPD:
    flux.ceair = min(max(flux.eair, 0.03*leaf.esat), leaf.esat )

    # Electron transport rate for C3 plants
    # Actual colimited potential Je (curvature and Jmax)
    flux.Je = minimum(quadratic(leaf.θ_j, -(flux.Je_red + leaf.Jmax), flux.Je_red * leaf.Jmax))    # Bonan eq. 11.21

    # Ci calculation
    # Medlyn or Ball-Berry:
    if leaf.dynamic_state # Save actual gs
        gs_actual = leaf.gs
    end

    if (leaf.gstyp <= 1)
        Ci_0 = leaf.C3 ? 0.7*met.Ca : 0.4*met.Ca
        # Solve iterative loop:
        leaf.Ci = hybrid(flux,leaf, met, CiFunc!, Ci_0, 1.05*Ci_0, tol)
    elseif leaf.gstyp == 2 # Needed for Bonan Stomatal optimization model
        leaf.Ci = CiFuncGs!(leaf.gs, flux, leaf, met)
    end
    if leaf.dynamic_state
        leaf.gs_ss = leaf.gs
        leaf.gs = gs_actual
        leaf.Ci = CiFuncGs!(leaf.gs, flux, leaf, met)
    end

    # Rate of actual CO2 per electron, incl. photorespiration
    # (Ci-Gamma_star)./(Ci+2*Gamma_star)
    leaf.CO2_per_electron = (leaf.Ci-leaf.Γstar)/(leaf.Ci+2leaf.Γstar) * leaf.effcon;

    # Actual effective ETR:
    flux.Ja = max(0,flux.Ag / leaf.CO2_per_electron);
    flux.Ja = min(flux.Ja,flux.Je_pot )

    # Effective photochemical yield:
    flux.φ = leaf.maxPSII*flux.Ja/flux.Je_pot;
    #println(flux.Ja, " ", flux.Je_pot)
    flux.φ = min(1/leaf.maxPSII,flux.φ)
    x   = max(0,  1-flux.φ/leaf.maxPSII);       # degree of light saturation: 'x' (van der Tol e.Ap. 2014)
    Fluorescencemodel!(flux.φ,x,leaf)


end # LeafPhotosynthesis (similar to biochem in SCOPE)









"""
    CiFunc!(Ci::Number, flux::fluxes, leaf::leaf_params, met::meteo)

Compute Assimilation using Ci as input

# Arguments
- `Ci::Number`: Ci.
- `flux::fluxes`: fluxes structure.
- `leaf::leaf_params`: leaf_params structure.
"""
function CiFunc!(Ci::Number, flux::fluxes, leaf::leaf_params, met::meteo)

    if leaf.C3
        # C3: Rubisco-limited photosynthesis; still need to check CO2 mixing ratios vs partial pressures.
        # still need to include ppm2bar (but can be done on leaf structure!)
        flux.Ac = leaf.Vcmax * max(Ci-leaf.Γstar, 0.0) / (Ci + leaf.Kc*(1.0+leaf.o₂/leaf.Ko)) # Bonan eq. 11.28
        # C3: RuBP-limited photosynthesis (this is the NADPH requirement stochiometry)
        flux.Aj = flux.Je * max(Ci-leaf.Γstar, 0.0) / (4.0*Ci + 8.0*leaf.Γstar)               # Bonan eq. 11.29

        # for C3, set ap to Inf
        flux.Ap = Inf
    else #C4 Photosynthesis, still to be implemented
        flux.Ac = flux.Aj = flux.Ap = 0.0
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
    flux.Ag = max(0.0,flux.Ag)
    flux.Ai = max(0.0,flux.Ai)
    flux.Aj = max(0.0,flux.Aj)
    flux.Ap = max(0.0,flux.Ap)

    # Net photosynthesis
    flux.An = flux.Ag - leaf.Rdleaf # net assimilation

    # adjust aerodynamic resistance based on leaf boundary layer and Monin Obukhov
    setra!(leaf, flux, met)

    # CO2 at leaf surface # might need to be changed
    flux.Cs = met.Ca - flux.An * leaf.ra/flux.g_m_s_to_micromol_m2_s

                                                          """
                                                          # Pierre to Christian: I really don't understand this part - I commented
                                                          # Stomatal constraint function (not sure we "need" the quadratic colimitations here, why not just use BB or Medlyn?)
                                                          if (leaf.gstyp == 1) # Ball-Berry
                                                              if flux.An >0.0
                                                                  leaf.gs = maximum(quadratic(flux.Cs, flux.Cs*(flux.gbv - leaf.g0) - leaf.g1*flux.An, -flux.gbv * (flux.Cs*leaf.g0 + leaf.g1*flux.An*flux.ceair/leaf.esat)))
                                                                  leaf.g1 * flux.An * flux.ceair/leaf.esat/flux.Cs  + leaf.g0;
                                                                  # println(leaf.gs)
                                                              else
                                                                  leaf.gs = leaf.g0
                                                              end
                                                          elseif (leaf.gstyp == 0) # Medlyn
                                                              if flux.An >0.0
                                                                  # Not sure how this all works, copied from Bonan's ML canopy model
                                                                  vpd_term = max((leaf.esat - flux.ceair), vpd_min) * 0.001
                                                                  term = 1.6 * flux.An / flux.Cs
                                                                  leaf.gs = maximum(quadratic(1.0, -(2.0 * (leaf.g0 + term) + (leaf.g1 * term)^2 / (flux.gbv * vpd_term)), leaf.g0 * leaf.g0 + (2.0 * leaf.g0 + term * (1.0 - leaf.g1 * leaf.g1 / vpd_term)) * term))
                                                              else
                                                                  leaf.gs = leaf.g0
                                                              end
                                                          end
                                                          """



    # compute stomatal conductance gs
    leaf.VPD       = max(leaf.esat-met.e_air,1.0); # can be negative at spin up
    leaf.RH        = min(max(met.e_air/leaf.esat,0.001),0.999);    # will need to be corrected alter to define surface RH

    if (leaf.gstyp == 1)
        BallBerry!(flux, leaf)
    else # Medlyn default
        Medlyn!(flux, leaf)
    end

    # Diffusion (supply-based) photosynthetic rate - Calculate Ci from the diffusion rate
    # total conductance
    leaf.gleaf = 1.0 / (leaf.ra/flux.g_m_s_to_micromol_m2_s + 1.6/leaf.gs + 1.0/leaf.gm)
    Cinew = met.Ca - flux.An / leaf.gleaf

    # CiFunc returns the difference between the current Ci and the new Ci
    leaf.Ci = Cinew
    return flux.An<0. ? 0.0 : Cinew - Ci
end




"""
    CiFuncGs!(gs::Number, flux::fluxes, leaf::leaf_params, met::meteo)

Compute Assimilation using fixed stomatal conductance gs.

# Arguments
- `gs::Number`: Stomatal conductance.
- `flux::fluxes`: fluxes structure.
- `leaf::leaf_params`: leaf_params structure.
"""
function CiFuncGs!(gs::Number, flux::fluxes, leaf::leaf_params, met::meteo)
    # adjust aerodynamic resistance based on leaf boundary layer and Monin Obukhov
    setra!(leaf, flux, met)

    # Compute overall conductance (Boundary layer, stomata and mesophyll)
    leaf.gleaf = 1.0/(leaf.ra/flux.g_m_s_to_micromol_m2_s + 1.6/gs + 1.0/leaf.gm)
    if leaf.gleaf<eps() leaf.gleaf=eps() end

    if leaf.C3
        # C3 Rubisco Limited Photosynthesis co-limited by gs
        a0 = leaf.Vcmax
        e0 = 1.0
        d0 = leaf.Kc*(1.0+leaf.o₂/leaf.Ko)
        flux.Ac = minimum(quadratic(1.0/leaf.gleaf, -(e0*met.Ca + d0) - (a0 - e0*leaf.Rdleaf) / leaf.gleaf, a0 * (met.Ca - leaf.Γstar) - leaf.Rdleaf * (e0*met.Ca + d0)))

        # C3: RuBP-limited photosynthesis
        a0 = flux.Je
        e0 = 4.0
        d0 = 8.0*leaf.Γstar
        flux.Aj = minimum(quadratic(e0 / leaf.gleaf, -(e0*met.Ca + d0) - (a0 - e0*leaf.Rdleaf) / leaf.gleaf, a0 * (met.Ca - leaf.Γstar) - leaf.Rdleaf * (e0*met.Ca + d0)))

        # C3: Product-limited photosynthesis
        flux.Ap = Inf
    # C4 to be implemented
    elseif !leaf.C3
        flux.Ac = flux.Aj = flux.Ap = 0.0
    end
    if leaf.use_colim
        flux.Ai = minimum(quadratic(leaf.C3 ? 0.98 : 0.80, -(flux.Ac + flux.Aj), flux.Ac * flux.Aj))   # Bonan Eq 11.33
        # Ap limitation only for C4 here:
        if leaf.C3
            flux.Ag = flux.Ai
        else # C4 colimitation with ap
            flux.Ag = minimum(quadratic(0.95, -(flux.Ai + flux.Ap), flux.Ai*flux.Ap))                  # Bonan Eq 11.33
        end
    else
        flux.Ag = min(flux.Ac,flux.Aj,flux.Ap)
    end
    flux.Ag = max(0,flux.Ag)
    flux.Ai = max(0,flux.Ai)
    flux.Aj = max(0,flux.Aj)
    flux.Ap = max(0,flux.Ap)

    # Compute net Photosynthesis
    flux.An = flux.Ag - leaf.Rdleaf

    # Compute CO2 at leaf surface
    flux.Cs = met.Ca - flux.An * leaf.ra/flux.g_m_s_to_micromol_m2_s

    # Compute Ci (included Mesophyll as well in principle)
    ci_val = met.Ca - flux.An /leaf.gleaf

    #leaf.CO2_per_electron = (ci_val-leaf.Γstar)./(ci_val+2.0*leaf.Γstar) .* leaf.effcon;
end # Function CiFuncGs!



function ψ_m(ζ) # momentum correction function
    # stability corrections - Zeng et al., 1998
    #if(ζ>1.0) # very stable - Grachev et al SHEBA 2007

    if(ζ>=0.0) # stable - Grachev et al SHEBA 2007
        x  = (1.0+ζ)^(1.0/3.0);
        am = 5.0;
        bm = am/6.5;
        ah = bh = 5.0;
        ch = 3.0;
        Bm = ((1.0-bm)/bm)^(1.0/3.0);
        Bh = sqrt(5);
        ψ = -3.0*am/bm*(x-1.0) + am*Bm/(2.0*bm)*( 2.0*log((x+Bm)/(1.0+Bm))
                                                 - log( (x*x-x*Bm+Bm*Bm)/(1.0-Bm+Bm*Bm) )
                                                 +2*sqrt(3.0)*( atan((2.0*x-Bm)/(sqrt(3.0)*Bm)) -  atan((2.0-Bm)/(sqrt(3.0)*Bm)))  );
    #elseif (ζ<=1.0 and ζ>=0)  # stable
    #    ψ = -5.0*ζ;
    #    ψ_h = -5.0*ζ;

    elseif (ζ<0 && ζ>=-0.465) # unstable
        x   = (1.0-16.0*ζ)^0.25;
        ψ = 2.0*log((1.0+x)/2.0) + log((1.0+x*x)/2.0) - 2.0*atan(x) + π/2;
    elseif (ζ<-0.465) # very unstable
        x   = (1.0-16.0*ζ)^0.25;
        ψ = 2.0*log((1.0+x)/2.0) + log((1.0+x*x)/2.0) - 2.0*atan(x) + π/2;
    end
    return ψ
end

function ψ_h(ζ) # momentum correction function
    # stability corrections - Zeng et al., 1998
    #if(ζ>1.0) # very stable - Grachev et al SHEBA 2007

    if(ζ>=0.0) # stable - Grachev et al SHEBA 2007
        x  = (1.0+ζ)^(1.0/3.0);
        am = 5.0;
        bm = am/6.5;
        ah = bh = 5.0;
        ch = 3.0;
        Bm = ((1.0-bm)/bm)^(1.0/3.0);
        Bh = sqrt(5);
        ψ = -bh/2*log(1+ch*ζ+ζ*ζ) + (-ah/Bh + bh*ch/(2.0*Bh))*(log( (2.0*ζ+ch-Bh)/(2.0*ζ+ch+Bh) - log((ch-Bh)/(ch+Bh)) ));
    #elseif (ζ<=1.0 and ζ>=0)  # stable
    #    ψ_m = -5.0*ζ;
    #    ψ = -5.0*ζ;

    elseif (ζ<0 && ζ>=-0.465) # unstable
        x   = (1.0-16.0*ζ)^0.25;
        ψ = 2.0*log((1.0+x*x)/2.0);
    elseif (ζ<-0.465) # very unstable
        x   = (1.0-16.0*ζ)^0.25;
        ψ = 2.0*log((1.0+x*x)/2.0);
    end
    return ψ
end

function setra_atmo(l::leaf_params,flux::fluxes,met::meteo)
    # based on Monin-Obukhov Similiarity theory -> to be changed for LES
    # compute Obukhov length
    Hv_s  = flux.H*(1+0.61*met.e_air);
    Tv    = met.T_air*(1+0.61*met.e_air);
    L     = - flux.ustar^3*Tv/(physcon.grav*physcon.K*Hv_s);
    ra_m = 1.0/(physcon.K^2*met.U) * ( log((met.zscreen - l.d)/l.z0m) - ψ_m((met.zscreen - l.d)/L) + ψ_m(l.z0m/L) ) * ( log((met.zscreen - l.d)/l.z0m) - ψ_m((met.zscreen - l.d)/L) + ψ_h(l.z0h/L) ) ;# momentum aerodynamic resistance
    ra_w = 1.0/(physcon.K^2*met.U) * ( log((met.zscreen - l.d)/l.z0m) - ψ_m((met.zscreen - l.d)/L) + ψ_m(l.z0m/L) ) * ( log((met.zscreen - l.d)/l.z0h) - ψ_h((met.zscreen - l.d)/L) + ψ_h(l.z0h/L) ) ;# water aerodynamic resistance
    flux.ustar = sqrt(met.U/ra_m); # keep ustar estimate in memory
    return ra_w; # based on Monin-Obukhov Similarity theory -> to be changed for LES
end

function setra!(l::leaf_params, flux::fluxes, met::meteo) # set leaf boundary layer
    ra_leaf = 1/l.Cd / sqrt(met.U/l.dleaf); # kmax . int_psis^psil k(x)dx = kmax . IntWeibull(psil);
    l.ra    = ra_leaf + setra_atmo(l, flux, met);
end















"""
    Fluorescencemodel!(ps,x,leaf::leaf_params )

Compute Fluorescence yields, Kn and Kp.

# Arguments
- `ps::Number`: PSII yield.
- `x::Number`: Degree of light saturation: [0-1] .
- `leaf::leaf_params`: leaf_params structure.
"""
function Fluorescencemodel!(ps::Number,x::Number,leaf::leaf_params )
    x_alpha = exp(log(x)*leaf.Knparams[2]); # this is the most expensive operation in this fn; doing it twice almost doubles the time spent here (MATLAB 2013b doesn't optimize the duplicate code)
    #println(x_alpha)
    leaf.Kn_ss = leaf.Knparams[1] * (1+leaf.Knparams[3])* x_alpha/(leaf.Knparams[3] + x_alpha);
    Kf = leaf.Kf
    Kn = leaf.Kn
    Kd = leaf.Kd
    leaf.Kp   = max(0,-ps*(Kf+Kd+Kn)/(ps-1));
    Kp = leaf.Kp


    leaf.Fo   = Kf/(Kf+4.0+Kd);
    leaf.Fo′  = Kf/(Kf+4.0+Kd+Kn);
    leaf.Fm   = Kf/(Kf   +Kd);
    leaf.Fm′  = Kf/(Kf   +Kd+Kn);
    leaf.ϕs   = leaf.Fm′*(1-ps);
    leaf.eta  = leaf.ϕs/leaf.Fo;
    leaf.qQ   = 1-(leaf.ϕs-leaf.Fo′)/(leaf.Fm-leaf.Fo′);
    leaf.qE   = 1-(leaf.Fm-leaf.Fo′)/(leaf.Fm′-leaf.Fo);

    leaf.NPQ  = Kn/(Kf+Kd);
end









# mathematical functions

function hybrid(flux::fluxes,leaf::leaf_params, met::meteo, func::Function, xa::Number, xb::Number, tol::Number)
    #
    # !DESCRIPTION:
    # Solve for the root of a function, given initial estimates xa and xb.
    # The root is updated until its accuracy is tol.
    #
    # !USES:
    #
    # !ARGUMENTS:
    #procedure (xfunc) :: func         ! Function to solve
    #real(r8), intent(in) :: xa, xb    ! Initial estimates of root
    #real(r8), intent(in) :: tol       ! Error tolerance
    #type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    itmax = 40
    x0 = x1 = f0 = f1 = 0.0
    x0 = xa
    f0 = func(x0,flux, leaf, met)
    if (f0 == 0.0) return  x0 end

    x1 = xb
    f1 = func(x1,flux, leaf, met)
    if (f1 == 0.0) return x1 end

    if (f1 < f0)
       minx = x1
       minf = f1
    else
       minx = x0
       minf = f0
    end

    # First use the secant method, and then use the brent method as a backup
    for iter = 0:itmax
       iter = iter + 1
       dx = -f1 * (x1 - x0) / (f1 - f0)
       x = x1 + dx
       if (abs(dx) < tol)
          x0 = x
          break
       end
       x0 = x1
       f0 = f1
       x1 = x
       f1 = func(x1,flux, leaf, met)
       if (f1 < minf)
          minx = x1
          minf = f1
       end

       # If a root zone is found, use the brent method for a robust backup strategy
       if (f1 * f0 < 0.0)
          x = zbrent(flux, leaf, met, func,x0, x1, xtol=tol)
          x0 = x
          break
       end

       # In case of failing to converge within itmax iterations stop at the minimum function
       if (iter > itmax)
          f1 = func(minx,flux, leaf. met)
          x0 = minx
          break
       end

    end

    return x0
end #function hybrid

function zbrent(flux::fluxes,leaf::leaf_params,met::meteo,f::Function, x0::Number, x1::Number, args::Tuple=();
               xtol::AbstractFloat=1e-7, ytol=2eps(Float64),
               maxiter::Integer=50)
    EPS = eps(Float64)
    y0 = f(x0,flux,leaf)
    y1 = f(x1,flux,leaf)
    if abs(y0) < abs(y1)
        # Swap lower and upper bounds.
        x0, x1 = x1, x0
        y0, y1 = y1, y0
    end
    x2 = x0
    y2 = y0
    x3 = x2
    bisection = true
    for _ in 1:maxiter
        # x-tolerance.
        if abs(x1-x0) < xtol
            return x1
        end

        # Use inverse quadratic interpolation if f(x0)!=f(x1)!=f(x2)
        # and linear interpolation (secant method) otherwise.
        if abs(y0-y2) > ytol && abs(y1-y2) > ytol
            x = x0*y1*y2/((y0-y1)*(y0-y2)) +
                x1*y0*y2/((y1-y0)*(y1-y2)) +
                x2*y0*y1/((y2-y0)*(y2-y1))
        else
            x = x1 - y1 * (x1-x0)/(y1-y0)
        end

        # Use bisection method if satisfies the conditions.
        delta = abs(2EPS*abs(x1))
        min1 = abs(x-x1)
        min2 = abs(x1-x2)
        min3 = abs(x2-x3)
        if (x < (3x0+x1)/4 && x > x1) ||
           (bisection && min1 >= min2/2) ||
           (!bisection && min1 >= min3/2) ||
           (bisection && min2 < delta) ||
           (!bisection && min3 < delta)
            x = (x0+x1)/2
            bisection = true
        else
            bisection = false
        end

        y = f(x,flux,leaf)
        # y-tolerance.
        if abs(y) < ytol
            return x
        end
        x3 = x2
        x2 = x1
        if sign(y0) != sign(y)
            x1 = x
            y1 = y
        else
            x0 = x
            y0 = y
        end
        if abs(y0) < abs(y1)
            # Swap lower and upper bounds.
            x0, x1 = x1, x0
            y0, y1 = y1, y0
        end
    end
    error("Max iteration exceeded")
end

include("../Utils/math_tools.jl")


end #Module
