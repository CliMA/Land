module LeafPhotosynthesisMod

using Parameters

using ..PhysCon
using ..WaterVaporMod

include("Leaf.jl")
include("surface_energy_balance.jl")

export fluxes, ψ_m, ψ_h, setra!
"Tolerance thhreshold for Ci iterations"
tol = 0.1
vpd_min = 0.1

" Just a placeholder for now"
@with_kw mutable struct fluxes{TT<:Number}
         APAR::TT = 500.0
         gbc::TT = 100.0
         gbv::TT = 100.0
         cair::TT = 400.0
         ceair::TT = 1400.0
         eair::TT = 1400.0
         je::TT = 1100.0
         #gs::TT = 0.0
         ac::TT = 0.0
         aj::TT = 0.0
         ai::TT = 0.0
         ap::TT = 0.0
         ag::TT = 0.0
         an::TT = 0.0
         cs::TT = 0.0
         #ci::TT = 0.0
         rd::TT = 0.0
         Je_pot::TT = 0.0
         Ja::TT = 0.0
         Je_red::TT = 0.0
         φ::TT = 0.0
         U::TT = 0.0
end

struct atmos
         o2air              # Atmospheric O2 (mmol/mol)
         co2air             # Atmospheric CO2 (μmol/mol)
         eair               # Vapor pressure of air (Pa)
end

"""
    LeafPhotosynthesis(flux::fluxes, leaf::leaf_params,T::Number)

Compute net assimilation rate A, fluorescence F using biochemical model

# Arguments
- `flux::fluxes`: fluxes structure.
- `leaf::leaf_params`: leaf_params structure.
- `T::Number`: Leaf Temperature
"""
function LeafPhotosynthesis(flux::fluxes, leaf::leaf_params,T::Number)
    # Adjust rates to leaf Temperature (C3 only for now):
    setLeafT!(leaf, T)
    # adjust aerodynamic resistance based on leaf boundary layer and Monin Obukhov
    setra!(leaf, flux.U )

    # Compute max PSII efficiency here (can later be used with a variable Kn!)
    leaf.Kp = 4.0
    φ_PSII  = leaf.Kp/(leaf.Kp+leaf.Kf+leaf.Kd+leaf.Kn)

    # Save leaf respiration
    flux.rd = leaf.rdleaf;
    # Calculate potential electron transport rate (assuming no upper bound, proportional to absorbed light!):
    flux.Je_pot = 0.5 * leaf.maxPSII * flux.APAR;                          # potential electron transport rate (important for later)
    flux.Je_red = 0.5 * φ_PSII * flux.APAR;                                # Includes Kn here
    # Some bound constraint on VPD:
    flux.ceair = min(max(flux.eair, 0.03*leaf.esat), leaf.esat )

    # Electron transport rate for C3 plants
    # Actual colimited potential Je (curvature and Jmax)
    flux.je = minimum(quadratic(leaf.θ_j, -(flux.Je_red + leaf.jmax), flux.Je_red * leaf.jmax))    # Bonan eq. 11.21

    # Ci calculation
    # Medlyn or Ball-Berry:
    if leaf.dynamic_state # Save actual gs
        gs_actual = leaf.gs
    end

    if (leaf.gstyp <= 1)
        Ci_0 = leaf.C3 ? 0.7*flux.cair : 0.4*flux.cair
        # Solve iterative loop:
        leaf.Ci = hybrid(flux,leaf, CiFunc!, Ci_0, 1.05*Ci_0, tol)
    elseif leaf.gstyp == 2 # Needed for Bonan Stomatal optimization model
        leaf.Ci = CiFuncGs!(leaf.gs, flux,leaf)
    end
    if leaf.dynamic_state
        leaf.gs_ss = leaf.gs
        leaf.gs = gs_actual
        leaf.Ci = CiFuncGs!(leaf.gs, flux,leaf)
    end

    # Rate of actual CO2 per electron, incl. photorespiration
    # (Ci-Gamma_star)./(Ci+2*Gamma_star)
    leaf.CO2_per_electron = (leaf.Ci-leaf.Γstar)/(leaf.Ci+2leaf.Γstar) * leaf.effcon;

    # Actual effective ETR:
    flux.Ja = max(0,flux.ag / leaf.CO2_per_electron);
    flux.Ja = min(flux.Ja,flux.Je_pot )

    # Effective photochemical yield:
    flux.φ = leaf.maxPSII*flux.Ja/flux.Je_pot;
    #println(flux.Ja, " ", flux.Je_pot)
    flux.φ = min(1/leaf.maxPSII,flux.φ)
    x   = max(0,  1-flux.φ/leaf.maxPSII);       # degree of light saturation: 'x' (van der Tol e.a. 2014)
    Fluorescencemodel!(flux.φ,x,leaf)


end # LeafPhotosynthesis (similar to biochem in SCOPE)









"""
    CiFunc!(Ci::Number, flux::fluxes, leaf::leaf_params)

Compute Assimilation using Ci as input

# Arguments
- `Ci::Number`: Ci.
- `flux::fluxes`: fluxes structure.
- `leaf::leaf_params`: leaf_params structure.
"""
function CiFunc!(Ci::Number, flux::fluxes, leaf::leaf_params)

    if leaf.C3
        # C3: Rubisco-limited photosynthesis; still need to check CO2 mixing ratios vs partial pressures.
        # still need to include ppm2bar (but can be done on leaf structure!)
        flux.ac = leaf.vcmax * max(Ci-leaf.Γstar, 0.0) / (Ci + leaf.kc*(1.0+leaf.o₂/leaf.ko)) # Bonan eq. 11.28
        # C3: RuBP-limited photosynthesis (this is the NADPH requirement stochiometry)
        flux.aj = flux.je * max(Ci-leaf.Γstar, 0.0) / (4.0*Ci + 8.0*leaf.Γstar)               # Bonan eq. 11.29

        # for C3, set ap to Inf
        flux.ap = Inf
    else #C4 Photosynthesis, still to be implemented
        flux.ac = flux.aj = flux.ap = 0.0
    end
    # Net photosynthesis as the minimum or co-limited rate
    if leaf.use_colim
        flux.ai = minimum(quadratic(leaf.C3 ? 0.98 : 0.80, -(flux.ac + flux.aj), flux.ac * flux.aj))
        if leaf.C3
            flux.ag = flux.ai
        else # C4 colimitation with ap
            flux.ag = minimum(quadratic(0.95, -(flux.ai + flux.ap), flux.ai*flux.ap))                 # Bonan Eq 11.33
        end
    else
        flux.ag = min(flux.ac,flux.aj,flux.ap)
    end
    # Prevent photosynthesis from ever being negative
    flux.ag = max(0,flux.ag)
    flux.ai = max(0,flux.ai)
    flux.aj = max(0,flux.aj)
    flux.ap = max(0,flux.ap)

    # Net photosynthesis
    flux.an = flux.ag - leaf.rdleaf

    # CO2 at leaf surface # might need to be changed
    flux.cs = flux.cair - flux.an / flux.gbc

    # Stomatal constraint function (not sure we "need" the quadratic colimitations here, why not just use BB or Medlyn?)
    if (leaf.gstyp == 1) # Ball-Berry
        if flux.an >0.0
            leaf.gs = maximum(quadratic(flux.cs, flux.cs*(flux.gbv - leaf.g0) - leaf.g1*flux.an, -flux.gbv * (flux.cs*leaf.g0 + leaf.g1*flux.an*flux.ceair/leaf.esat)))
            leaf.g1 * flux.an * flux.ceair/leaf.esat/flux.cs  + leaf.g0;
            # println(leaf.gs)
        else
            leaf.gs = leaf.g0
        end
    elseif (leaf.gstyp == 0) # Medlyn
        if flux.an >0.0
            # Not sure how this all works, copied from Bonan's ML canopy model
            vpd_term = max((leaf.esat - flux.ceair), vpd_min) * 0.001
            term = 1.6 * flux.an / flux.cs
            leaf.gs = maximum(quadratic(1.0, -(2.0 * (leaf.g0 + term) + (leaf.g1 * term)^2 / (flux.gbv * vpd_term)), leaf.g0 * leaf.g0 + (2.0 * leaf.g0 + term * (1.0 - leaf.g1 * leaf.g1 / vpd_term)) * term))
        else
            leaf.gs = leaf.g0
        end
    end
    # Diffusion (supply-based) photosynthetic rate - Calculate Ci from the diffusion rate
    gleaf = 1.0 / (1.0/flux.gbc + 1.6/leaf.gs + 1.0/leaf.gm)
    cinew = flux.cair - flux.an / gleaf

    # CiFunc returns the difference between the current Ci and the new Ci
    leaf.Ci = cinew
    return flux.an<0. ? 0.0 : cinew - Ci
end




"""
    CiFuncGs!(gs::Number, flux::fluxes, leaf::leaf_params)

Compute Assimilation using fixed stomatal conductance gs.

# Arguments
- `gs::Number`: Stomatal conductance.
- `flux::fluxes`: fluxes structure.
- `leaf::leaf_params`: leaf_params structure.
"""
function CiFuncGs!(gs::Number, flux::fluxes, leaf::leaf_params)
    # Compute overall conductance (Boundary layer, stomata and mesophyll)
    gleaf = 1.0/(1.0/flux.gbc + 1.6/gs + 1.0/leaf.gm)
    if gleaf<eps() gleaf=eps() end

    if leaf.C3
        # C3 Rubisco Limited Photosynthesis co-limited by gs
        a0 = leaf.vcmax
        e0 = 1.0
        d0 = leaf.kc*(1.0+leaf.o₂/leaf.ko)
        flux.ac = minimum(quadratic(1.0/gleaf, -(e0*flux.cair + d0) - (a0 - e0*leaf.rdleaf) / gleaf, a0 * (flux.cair - leaf.Γstar) - leaf.rdleaf * (e0*flux.cair + d0)))

        # C3: RuBP-limited photosynthesis
        a0 = flux.je
        e0 = 4.0
        d0 = 8.0*leaf.Γstar
        flux.aj = minimum(quadratic(e0 / gleaf, -(e0*flux.cair + d0) - (a0 - e0*leaf.rdleaf) / gleaf, a0 * (flux.cair - leaf.Γstar) - leaf.rdleaf * (e0*flux.cair + d0)))

        # C3: Product-limited photosynthesis
        flux.ap = Inf
    # C4 to be implemented
    elseif !leaf.C3
        flux.ac = flux.aj = flux.ap = 0.0
    end
    if leaf.use_colim
        flux.ai = minimum(quadratic(leaf.C3 ? 0.98 : 0.80, -(flux.ac + flux.aj), flux.ac * flux.aj))   # Bonan Eq 11.33
        # Ap limitation only for C4 here:
        if leaf.C3
            flux.ag = flux.ai
        else # C4 colimitation with ap
            flux.ag = minimum(quadratic(0.95, -(flux.ai + flux.ap), flux.ai*flux.ap))                  # Bonan Eq 11.33
        end
    else
        flux.ag = min(flux.ac,flux.aj,flux.ap)
    end
    flux.ag = max(0,flux.ag)
    flux.ai = max(0,flux.ai)
    flux.aj = max(0,flux.aj)
    flux.ap = max(0,flux.ap)
    # Compute net Photosynthesis
    flux.an = flux.ag - leaf.rdleaf
    # Compute CO2 at leaf surface
    flux.cs = flux.cair - flux.an / flux.gbc

    # Compute Ci (included Mesophyll as well in principle)
    ci_val = flux.cair - flux.an / gleaf
    #leaf.CO2_per_electron = (ci_val-leaf.Γstar)./(ci_val+2.0*leaf.Γstar) .* leaf.effcon;
end # Function CiFuncGs!






function setra!(l::leaf_params, flux::fluxes, met::meteo) # set leaf boundary layer
    ra_leaf = 1/l.Cd / sqrt(met.U/l.dleaf); # kmax . int_psis^psil k(x)dx = kmax . IntWeibull(psil);
    l.ra    = ra_leaf + setra_atmo(l,flux, met);
end


function setra_atmo!(l::leaf_params,flux::fluxes,met::meteo)
    # based on Monin-Obukhov Similiarity theory -> to be changed for LES
    # compute Obukhov length
    flux.Hv_s = flux.H*(1+0.61*met.ea);
    L = - ustar^3*Tv/(physcon.grav*physcon.K*Hv_s);
    ra_w = 1.0/(physcon.K^2*met.U) * ( log((met.zscreen - l.d)/met.z0m) - ψ_m((met.zscreen - l.d)/L) + ψ_m(met.z0m/L) ) * ( log((met.zscreen - l.d)/met.z0h) - ψ_h((met.zscreen - l.d)/L) + ψ_h(met.z0h/L) ) ;# water aerodynamic resistance
    return ra_w; # based on Monin-Obukhov Similarity theory -> to be changed for LES
end

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

function hybrid(flux::fluxes,leaf::leaf_params, func::Function, xa::Number, xb::Number, tol::Number)
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
    f0 = func(x0,flux, leaf)
    if (f0 == 0.0) return  x0 end

    x1 = xb
    f1 = func(x1,flux, leaf)
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
       f1 = func(x1,flux, leaf)
       if (f1 < minf)
          minx = x1
          minf = f1
       end

       # If a root zone is found, use the brent method for a robust backup strategy
       if (f1 * f0 < 0.0)
          x = zbrent(flux,  leaf,func,x0, x1, xtol=tol)
          x0 = x
          break
       end

       # In case of failing to converge within itmax iterations stop at the minimum function
       if (iter > itmax)
          f1 = func(minx,flux, leaf)
          x0 = minx
          break
       end

    end

    return x0
end #function hybrid

function zbrent(flux::fluxes,leaf::leaf_params,f::Function, x0::Number, x1::Number, args::Tuple=();
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
