abstract type AbstractLeafBoundaryLayer end

"""
Gentine Boundary Layer computation, compute Monin Obhukow length and resistances across the leaf to Ca scale
"""
struct GentineLeafBoundary  <: AbstractLeafBoundaryLayer end


"""
Fixed resistance leaf Boundary Layer computation, ensures surface concentration [`CanopyLayer`](@ref)
"""
Base.@kwdef struct FixedBoundaryResistance{FT} <: AbstractLeafBoundaryLayer
  ra::FT = 0.0;
end

"""
    boundary_layer_resistance!(mod::NoResistanceBoundary, leaf, meteo)

Apply fixed leaf boundary layer resistance, given
- `mod` of type 
- `leaf` A [`CanopyLayer`](@ref) in [`Canopy`](@ref) in a [`Tree`](@ref)
- `meteo` indx'th [`Leaf`](@ref) in the `canopyi.leaf_list`
"""
function boundary_layer_resistance!(mod::FixedBoundaryResistance, l::leaf_params,  met::meteo) # set aerodynamic resistance
    # need to compute the evolving buoyancy flux
    Tv       = met.T_air*(1.0+0.61*met.e_air);

    #println("ra_leaf=",ra_leaf)
    ΔT   =   l.T - met.T_air;
    VPD      =   max(l.esat-met.e_air,1.0);
    #@show VPD
    ρd       =   met.P_air/(physcon.Rd*met.T_air);      # dry air density (kg/m3)
    lv       =   Lv(l.T);
    raw_full =   mod.ra;
    l.H      =   ρd*physcon.Cpd*ΔT/raw_full;
    rs_s_m   =   met.g_m_s_to_mol_m2_s/l.gs;
    l.LE     =   physcon.ε/met.P_air*ρd*lv*VPD/(rs_s_m+raw_full);
    met.ra = raw_full
end

"""
    get_marginal_gain(canopyi, indx)

Marginal water use efficiency `∂A/∂E` for a leaf in a canopy layer, given
- `canopyi` A [`CanopyLayer`](@ref) in [`Canopy`](@ref) in a [`Tree`](@ref)
- `indx` indx'th [`Leaf`](@ref) in the `canopyi.leaf_list`
"""
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