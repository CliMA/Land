
module WaterVapor

using ..PhysCon

#-----------------------------------------------------------------------
# DESCRIPTION:
# Calculate saturation vapor pressure and latent heat of vaporization
#
export Lv, SatVap, μ_l


function Lv(T)
  # latent heat of vaporization as a function of temperature
  # B. HENDERSON-SELLERS (1984) QJRMS
  # A new formula for latent heat of vaporization of water as a function of temperature
  Lv = 1.91846e6*(T/(T - 33.91))^2 #(J/kg)
  return Lv
  #es(T) = 2.1718e10 * exp(-4157/( T - 33.91))
end # end latent heat of vaporization

function SatVap(t)
    #
    # !DESCRIPTION:
    # Compute saturation vapor pressure and change in saturation vapor pressure
    # with respect to temperature. Polynomial approximations are from:
    # Flatau et al (1992) Polynomial fits to saturation vapor pressure.
    # Journal of Applied Meteorology 31:1507-1513
    #
    # !USES:
    # tfrz = 273.15;                # Freezing point of water (K)
    #
    # !ARGUMENTS:

    # t        ! Temperature (K)
    # es       ! Vapor pressure (Pa)
    # d(es)/d(t) (Pa/K)

    # For water vapor (temperature range is 0C to 100C)
    a0 =  6.11213476;
    a1 =  0.444007856;
    a2 =  0.143064234e-01;
    a3 =  0.264461437e-03;
    a4 =  0.305903558e-05;
    a5 =  0.196237241e-07;
    a6 =  0.892344772e-10;
    a7 = -0.373208410e-12;
    a8 =  0.209339997e-15;

    # and for derivative
    b0 =  0.444017302;
    b1 =  0.286064092e-01;
    b2 =  0.794683137e-03;
    b3 =  0.121211669e-04;
    b4 =  0.103354611e-06;
    b5 =  0.404125005e-09;
    b6 = -0.788037859e-12;
    b7 = -0.114596802e-13;
    b8 =  0.381294516e-16;

    # For ice (temperature range is -75C to 0C)
    c0 =  6.11123516;
    c1 =  0.503109514;
    c2 =  0.188369801e-01;
    c3 =  0.420547422e-03;
    c4 =  0.614396778e-05;
    c5 =  0.602780717e-07;
    c6 =  0.387940929e-09;
    c7 =  0.149436277e-11;
    c8 =  0.262655803e-14;

    # and for derivative
    d0 =  0.503277922;
    d1 =  0.377289173e-01;
    d2 =  0.126801703e-02;
    d3 =  0.249468427e-04;
    d4 =  0.313703411e-06;
    d5 =  0.257180651e-08;
    d6 =  0.133268878e-10;
    d7 =  0.394116744e-13;
    d8 =  0.498070196e-16;

    tc = t - physcon.tfrz
    # Not sure these are needed but why not (shouldn't really exceeed them)
    if (tc > 100.0) tc = 100.0 end
    if (tc < -75.0) tc = -75.0 end

    if (tc >= 0.0)
       es    = a0 + tc*(a1 + tc*(a2 + tc*(a3 + tc*(a4 + tc*(a5 + tc*(a6 + tc*(a7 + tc*a8)))))))
       desdt = b0 + tc*(b1 + tc*(b2 + tc*(b3 + tc*(b4 + tc*(b5 + tc*(b6 + tc*(b7 + tc*b8)))))))
    else
       es    = c0 + tc*(c1 + tc*(c2 + tc*(c3 + tc*(c4 + tc*(c5 + tc*(c6 + tc*(c7 + tc*c8)))))))
       desdt = d0 + tc*(d1 + tc*(d2 + tc*(d3 + tc*(d4 + tc*(d5 + tc*(d6 + tc*(d7 + tc*d8)))))))
    end

    es    = es    * 100.;            # Convert from mb to Pa
    desdt = desdt * 100.;            # Convert from mb to Pa
    return es, desdt

end # SatVap


#-----------------------------------------------------------------------
# DESCRIPTION:
# Liquid water depenendence on viscosity
# Reid, Prausnitz, & Poling (1987)

function μ_l(T)
  A   =   1.856e-14; 	     # Pa·s
  B   =   4209;            # K
  C   =  	0.04527;         # K-1
  D   = 	-3.376e-5;      # K-2
  μl  =   A*exp( B/T + C*T + D*T*T );
  return μl
end # end mu l



end # Module
