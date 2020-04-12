####
#### Surface energy balance
####

using LeafPhotosynthesisMod

export LeafEnergyWaterBalance

"""
# compute surface energy balance of a leaf Rn-H-LE=0
# Arguments
- `met::meteo`: meteorological forcing.
- `Ts_t::Float32`: Leaf temperature forcing.
- `T::Number`: Leaf Temperature
"""
function LeafEnergyWaterBalance(Tleaf, psileaf, met::meteo, l::leaf_params,  flux::fluxes, psi_s)
    # plugged Tleaf and psileaf outside of l as they are being iterated during fine time stepping - rest assumed constant
    # find solution to the surface energy LeafEnergyBalance, with a small thermal inertia: C dT/dt = Rn-H-LE and water balance of a leaf d.dH20/dt = Sap - T
    # first set parameter values that are assumed not to be changing with energy budget
    l.T            = Tleaf;
    l.psi_l        = psileaf;

    LeafPhotosynthesis(flux, l, met);

    #print(l.gs)

    flux.Rn      =   (1-l.α)*met.S_down +  met.L_down - l.ε*physcon.σ*Tleaf^4;
    #dRn_dTs = - 4*l.ε*physcon.σ*Tleaf^3;
    ρd            =  met.P_air/(physcon.Rd*met.T_air);      # dry air density (kg/m3)
    #println("ρd= ",ρd, "Cp= ",physcon.Cpd, "Tleaf=", Tleaf, " (K), Tair=",met.T_air, " (K)")
    flux.H        =  ρd*physcon.Cpd*(Tleaf - met.T_air)/l.ra;
    #dH_dT   = ρd*physcon.Cpd/l.ra;
    lv           =   Lv(Tleaf);
    rs_s_m       =   flux.g_m_s_to_micromol_m2_s/l.gs;
    flux.LE      =   physcon.ε/met.P_air*ρd*lv*l.VPD/(rs_s_m+l.ra);
    #dLE_dT      =   physcon.ε/met.P_air*ρd*lv*desat_dT/(rs_s_m+l.ra);
    setkx!(l,psi_s, psileaf) ;# set hydraulic conductivity as a function of psis and psi_l
    flux.Sap     =   (psi_s - l.psi_l)*l.kx # equal to int k(psi)dpsi;
    dT_dt        =   (flux.Rn-flux.H-flux.LE)/(l.LMA*l.c_leaf);
    dH2Ol_dt     =   (flux.Sap-flux.LE/lv)/l.Ctree;

    #print(flux.Cs, " ppm, ". l.VPD/1000.0, " (kPA), ", flux.An, " micromol/s/m2  ")
    println("Sdown= " , met.S_down, "W/m2, Rn=",flux.Rn,"W/m2, SEB=",flux.Rn-flux.H-flux.LE,"W/m2, H= ",flux.H, "W /m2, LE= ",flux.LE, "W /m2, dT_dt=",dT_dt*3600," (K/hr), rs=",rs_s_m, " (s/m) , ra=",l.ra, " (s/m) ")
    return dT_dt,dH2Ol_dt , flux.Rn, flux.H, flux.LE
end
