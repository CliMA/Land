####
#### Surface energy balance
####

export meteo, LeafEnergyWaterBalance

@with_kw mutable struct meteo{TT<:Number}
     S_down::TT = -999.;
     L_down::TT = -999.;
     T_air::TT  = -999.;
     ea_air::TT = -999.;
     P_air::TT  =  1e5 ;      # surface pressure (Pa)
     #PAR:TT     = -999.;
     Ca::TT     =  400.;
     PAR::TT    = -999.;
     U::TT      = 1e-6;
end

"""
# compute surface energy balance of a leaf Rn-H-LE=0
# Arguments
- `met::meteo`: meteorological forcing.
- `Ts_t::Float32`: Leaf temperature forcing.
- `T::Number`: Leaf Temperature
"""
function LeafEnergyWaterBalance(Tleaf, psileaf, met::meteo, l::leaf_params,  psi_s)
    # plugged Tleaf and psileaf outside of l as they are being iterated during fine time stepping - rest assumed constant
    # find solution to the surface energy LeafEnergyBalance, with a small thermal inertia: C dT/dt = Rn-H-LE and water balance of a leaf d.dH20/dt = Sap - T
    # first set parameter values that are assumed not to be changing with energy budget
    flux           = fluxes();
    flux.APAR      = met.PAR;
    flux.U         = met.U;
    l.T            = Tleaf;
    l.psi_l        = psileaf;
    #LeafPhotosynthesis(flux, l);

    A              = met.PAR;
    Cs             = met.Ca;
    # !!!!! temporary
    VPD           = max(l.esat-met.ea_air,1.0); # can be negative at spin up
    #print(Cs, " ppm, ". VPD, " (PA)), ", A, " micoml/s/m2,-  ", l)
    Medlyn!(Cs, VPD, A, l); # adjust conductance
    #print(l.gs)
    setkx!(l,psi_s, psileaf) ;# set hydraulic conductivity as a function of psis and psi_l

    Rn      =   (1-l.α)*met.S_down +  met.L_down - l.ε*physcon.σ*Tleaf^4;
    #dRn_dTs = - 4*l.ε*physcon.σ*Tleaf^3;
    ρd =  met.P_air/(physcon.Rd*met.T_air);      # dry air density (kg/m3)
    H  =  ρd*physcon.Cpd*(Tleaf - met.T_air)/l.ra;
    #dH_dT   = ρd*physcon.Cpd/l.ra;
    lv  =  Lv(Tleaf);
    LE          =   physcon.ε/met.P_air*ρd*lv*VPD/(1/l.gs+l.ra);
    #dLE_dT      =   physcon.ε/met.P_air*ρd*lv*desat_dT/(1/l.gs+l.ra);
    Sap         =   (psi_s - l.psi_l)*l.kx # equl toint k(psi)dpsi;
    dT_dt       =   Rn/(l.LMA*l.c_leaf);  #(Rn-H-LE)/(l.LMA*l.c_leaf);
    dH2Ol_dt    =   (Sap-LE/lv)/l.Ctree;
    #println("Rn=",Rn,"W/m2, SEB=",Rn-H-LE,"W/m2, H= ",H, "W /m2, LE= ",LE, "W /m2, dT_dt=",dT_dt*3600," (K/hr), gs=",l.gs,",ra=",l.ra,", ea=",met.ea_air,", esat=",esat,", dH2Ol_dt=",dH2Ol_dt)
    return dT_dt,dH2Ol_dt
end
