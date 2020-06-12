###############################################################################
#
# Grid information struct
# This struct should go to types.jl
#
###############################################################################
"""
    struct SPACGrid{FT,lat,lom,ele}

Information for a soil-plant-air continuum in a grid
"""
Base.@kwdef mutable struct SPACGrid{FT}
    # TODO generalize this function?
    # Grid information
    "The latitutional size `[degree]`"
    dim_lat::FT = FT(0.25)
    "The longitutional size `[degree]`"
    dim_lon::FT = FT(0.25)
    "The latitude `[degree]`"
    lat    ::FT = FT(40.0)
    "The longitude `[degree]`"
    lon    ::FT = FT(-100.0)
    "The elevation `[m]`"
    ele    ::FT = FT(0.0)
    "Atmospheric pressure `[Pa]`"
    p_atm  ::FT = FT(101325.0)
    "Atmospheric O₂ partial pressure"
    p_O₂   ::FT = FT(21176.925)

    # Soil and Air profile from the soil to top of canopy
    "Total number of vertical layers above ground"
    n_above::Int = 20
    "Total number of vertical layers below ground"
    n_below::Int = 10
    "A list of upper z for the veritical layers above ground `[m]`"
    zl_above    ::Array{FT,1} = collect(FT,  1.0: 1.0:20.0)
    "A list of lower z for the veritical layers below ground `[m]`"
    al_below    ::Array{FT,1} = collect(FT, -0.2:-0.2:-2.0)
    "An array of CO₂ partial pressure along the vertical layers `[Pa]`"
    p_CO₂_array ::Array{FT,1} = zeros(FT, n_above) + 40
    "An array of H₂O partial pressure along the vertical layers `[Pa]`"
    p_H₂O_array ::Array{FT,1} = zeros(FT, n_above) + 1500
    "An array of soil water potential `[MPa]`"
    p_soil_array::Array{FT,1} = zeros(FT, n_above) + FT(-0.5)
    "An array of soil water content"
    swc_array   ::Array{FT,1} = zeros(FT, n_above) + FT(0.4)
    "An array of air temperature `[K]`"
    t_air_array ::Array{FT,1} = zeros(FT, n_above) + FT(298.15)
    "An array of soil temperature `[K]`"
    t_soil_array::Array{FT,1} = zeros(FT, n_above) + FT(298.15)

    # Plant information
    "An array of trees in the grid"
    tree_array::Array = []
    "An array of tree density in the grid"
    density_array::Array{FT,1} = []
end





# TODO add a function to map plants into a SPACGrid struct,
# automatically break the tree canopy into afew layers as defined in the SPACGrid
function generate_canopy_and_root()
    ;
end
















"""
update_canopy_from_rt_module!(tree::Tree, canopy_rt::struct_canopy, canOpt_rt::struct_canopyOptProps, canRad_rt::struct_canopyRadiation)

Updates canopy information from the RT module, given
- `tree` A [`Tree`](@ref) type
- `canopy_rt` A [`struct_canopy`] type from the Plant module
- `cadOpt_rt` A [`struct_canopyOptProps`] type from the CanopyRT module
- `canRad_rt` A [`struct_canopyRadiation`] type from the CanopyRT module

This interface function is pending for leaf temperature...
"""
function update_canopy_from_rt_module!(tree::Tree, canopy_rt::Canopy4RT, canOpt_rt::AbstractCanopyOpti, canRad_rt::CanopyRadiation)
    # fraction of sunlit leaves in each layer
    fraction_sl = repeat(canopy_rt.lidf, outer=[ length(canopy_rt.lazitab) ]) / length(canopy_rt.lazitab)

    # update the PAR from canopyRT module to the Plant module
    if tree.n_canopy == canopy_rt.nlayers
        nlayers = canopy_rt.nlayers
        for pl_layer in 1:nlayers
            canopyi  = tree.canopy_list[pl_layer]
            rt_layer = nlayers + 1 - pl_layer
            
            # set the diffuse par to all the leaves
            canopyi.par_list .= canRad_rt.absPAR_shadeCab[rt_layer] * 1E6
            # add the direct PAR to sunlit leaves
            canopyi.par_list[1:end-1] .+= reshape(canRad_rt.absPAR_sunCab[:,:,rt_layer],(:,1))[:,1] * 1E6

            # calculate the fraction of sunlit and shaded leaves
            f_view = (canOpt_rt.Ps[rt_layer]+canOpt_rt.Ps[rt_layer+1]) / 2
            la_new = canopyi.la .* [f_view .* fraction_sl; 1-f_view]
            canopyi.f_view  = f_view
            canopyi.la_list = la_new

            # update leaf temperature
            canopyi.t_list[1:end-1] .= reshape(canRad_rt.T_sun3D[:,:,rt_layer],(:,1))[:,1]
            canopyi.t_list[end]      = canRad_rt.T_shade[rt_layer]

            # update leaf-to-air VPD
            canopyi.d_list = saturation_vapor_pressure.(canopyi.t_list) .- canopyi.p_H₂O
        end
    else
        println("Error: the canopy layer in the Tree differs from that in the struct_canopy!")
    end
end





#=
"""
# compute surface energy balance of a leaf Rn-H-LE=0
# Arguments
- `met::meteo`: meteorological forcing.
- `Ts_t::Float32`: Leaf temperature forcing.
- `T::Number`: Leaf Temperature
"""
#function LeafEnergyWaterBalance(mod::AbstractPhotosynthesis, Tleaf, psileaf, Cc, met::meteo, l::LeafPhotosynthesisParams,   psi_s)
function LeafEnergyWaterBalance(mod::AbstractPhotosynthesis, Tleaf, psileaf, Cc, met::MeteoParams, l::LeafParams,   psi_s)
    # plugged Tleaf and psileaf outside of l as they are being iterated during fine time stepping - rest assumed constant
    # find solution to the surface energy LeafEnergyBalance, with a small thermal inertia: C dT/dt = Rn-H-LE and water balance of a leaf d.dH20/dt = Sap - T
    # first set parameter values that are assumed not to be changing with energy budget
    l.T            = Tleaf;
    l.psi_l        = psileaf;
    l.Cc           = Cc;

    if(Tleaf<200.0 || psileaf>0.0)
        println("Error in Tleaf and psi_leaf - unphysical values")
        @show Tleaf
        @show psileaf
        process.exit(10)
    end
    LeafPhotosynthesis!(mod, l, met); # compute H, LE in there as well
    LeafPhotosynthesis!(mod, l, met);
    # Leaf response time in seconds (15min)
    tau_gs = 15.0*60.0
    # TODO aerodynamic resistance should be double for the sensible heat flux at leaf level, but do we really care?
    l.Cleaf = l.LMA* ( (1.0-l.RWC)*physcon.Cdryleaf + l.RWC* CP_L )/(1.0-l.RWC); # leaf conductivity
    # println("Cleaf=",l.Cleaf)
    lv           =   latent_heat_vapor(Tleaf);
    l.Rn      =   (1-l.α)*met.S_down +  2 * (met.L_down - l.ε* K_BOLTZMANN *Tleaf^4); # 2 is for two sides of the leaves
    #println("S_down=",met.S_down," , Ldown=",met.L_down," Tleaf=",Tleaf)
    #dRn_dTs = - 4*l.ε*physcon.σ*Tleaf^3;
    setkx!(l,psi_s, psileaf) ;# set hydraulic conductivity as a function of psis and psi_l
    l.kx = l.kx/1000
    #@show l.kx
    # TODO * viscosty or / viscosity?
    l.Sap     =   (psi_s - l.psi_l - ρ_H₂O * GRAVITY * l.height)*l.kx * relative_viscosity(Tleaf); # equal to int k(psi)dpsi + gravity term=-rho.g.mean(k)*height, includes also water viscosty;
    #println("Sap=",lv*l.Sap," W/m^2, kx=",l.kx,", rho.g.h=",physcon.ρw*physcon.grav*l.height)
    #@show flux.LE
    dT_dt        =   (l.Rn-l.H-l.LE)/l.Cleaf; # 2 times for up and down part of the leaves - TODO need to check this I am not sure I agree when integrated over the canopy
    dH2Ol_dt     =   (l.Sap-l.LE/lv)/l.Ctree*1000;
    dt           =   1.0; # one second time step for Cc - a bit arbitrary and not cruCcal but will be cahnged later for actual airspace
    #dCc_dt       =   0.0#(flux.An_diffusion - flux.An_biochemistry) / dt; #(l.Chloroplast_rel_volume*l.dleaf*l.LAI);
    dgs_dt       =   (l.gs_ss-l.gs)/tau_gs
    #println("Cc/(dCc/dt)=",1.0/(dCc_dt/l.Cc)," An_biochemistry=",flux.An_biochemistry," An_diffusion=",flux.An_diffusion)

    # print(flux.Cs, " ppm, ", l.VPD/1000.0, " (kPA), ", flux.An_biochemistry, " micromol/s/m2  "," An_diffusion=",flux.An_diffusion)
    #println("Sdown= " , met.S_down, "W/m2, Rn=",flux.Rn,"W/m2, SEB=",flux.Rn-flux.H-flux.LE,"W/m2, H= ",flux.H, "W /m2, LE= ",flux.LE, "W /m2, dT_dt=",dT_dt*3600," (K/hr), ra=",flux.ra, " (s/m) ")
    return dT_dt, dH2Ol_dt, dgs_dt
end
=#





#=
"""
  set_leaf_temperature!(mods::PhotoMods, l::LeafPhotosynthesisParams)

Computes Vcmax, Jmax, Rd, Kc, Ko, esat, desat, at given leaf temperature (if temperature changed), given
- `model` One [`PhotoMods`](@ref) model setup structure
- `leaf` One [`LeafPhotosynthesisParams`](@ref) structure 
"""
function set_leaf_temperature!(mod::AbstractPhotosynthesis, l::LeafPhotosynthesisParams)
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
=#




#=


# Go to SPAC later!

"""
    CcFunc!(mods::PhotoMods,  leaf::LeafPhotosynthesisParams, met::meteo)

Compute gross and net assimilation rates Ag, An using biochemical model and given Cc, using
- `model` One [`PhotoMods`](@ref) model setup structure
- `leaf` One [`LeafPhotosynthesisParams`](@ref) structure 
- `met` One [`lmeteo`](@ref) structure 

"""
function CcFunc!(Cc; mods::AbstractPhotosynthesis,  leaf::LeafPhotosynthesisParams, met::meteo)
    leaf.Cc = max(0,Cc)
    #@show leaf.Cc
    # Compute Rubisco Limited Photosynthesis
    rubisco_limited_rate!(mods.photosynthesis, leaf, met)
    # Light limited rate
    light_limited_rate!(mods.photosynthesis, leaf, met, leaf.APAR)
    # Product limited rate
    product_limited_rate!(mods.photosynthesis, leaf, met)
    
    # adjust aerodynamic resistance based on leaf boundary layer and Monin Obukhov
    boundary_layer_resistance!(mods.BoundaryLayer, leaf,  met)

    @unpack Ac, Aj,Ap, Cs, VPD, RH, gm, gs, Ag, Rd, An, gleaf, esat, Cc = leaf
    #@show Ac, Aj,Ap, Cs, VPD, RH, gm, gs, Ag, Rd, An, gleaf, esat, Cc 
    @unpack Ca, ra, g_m_s_to_mol_m2_s, e_air = met

    # add colimitation later:
    #@show Ac,Aj,Ap, Cc
    leaf.Ag = photosynthesis_colimit!(mods.colimitation, Ac, Aj, Ap) 

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


=#



#=



# Go to SPAC later!
"""
    LeafPhotosynthesis!(mod::AbstractPhotosynthesis, leaf::LeafPhotosynthesisParams, met::meteo)

Compute net assimilation rate A, fluorescence F using biochemical model, given
- `mod` One [`PhotoMods`](@ref) model setup structure
- `leaf` One [`LeafPhotosynthesisParams`](@ref) structure 
- `met` One [`lmeteo`](@ref) structure 

"""
function LeafPhotosynthesis!(mo::AbstractPhotosynthesis, leaf::LeafPhotosynthesisParams, met::meteo)
  # Set leaf temperature:
  set_leaf_temperature!(mo, leaf)
  met.g_m_s_to_mol_m2_s = (met.P_air-met.e_air)/(physcon.Rgas*met.T_air) 
  met.ppm_to_Pa = (met.P_air-met.e_air)*1e-6
  isnan(leaf.gleaf) ?  leaf.gleaf = 0.01 :
  isnan(leaf.gs) ?  leaf.gs = 0.01 :

  # Compute Cc and Photosynthesis:
  if leaf.dynamic_state
    CcFunc!(leaf.Cc; mods=mo, leaf=leaf, met=met)
  else
    #@show (CcFunc!(leaf.Cc; mods=mo, leaf=leaf, met=met))
    @inline f(x) = CcFunc!(x; mods=mo, leaf=leaf, met=met)
    sol = find_zero(f, SecantMethod{FT}(met.Ca*0.6, met.Ca), CompactSolution(),SolutionTolerance{FT}(1e-3))
    leaf.Cc = sol.root
    #@show sol
    #@show leaf.Cc
    #CcFunc!(mo, leaf, met)
    #@show leaf.Cc
  end

  # Compute Fluorescence
  leaf_fluorescence!(mo.fluorescence,  leaf);
  return leaf.An

end # LeafPhotosynthesis (similar to biochem in SCOPE)

=#






#=
"""
  CanopyPhotosynthesis!(mo::AbstractPhotosynthesis, leaf::LeafPhotosynthesisParams, met::meteo, canRad)

Compute Photosynthesis for all leaves in a canopy, given:
- `mod`    One [`PhotoMods`](@ref) model setup structure
- `leaf`   One [`LeafPhotosynthesisParams`](@ref) structure 
- `met`    One [`lmeteo`](@ref) structure 
- `canRad` One [`struct_canopyRadiation`](@ref) structure 
"""
function CanopyPhotosynthesis!(mo::AbstractPhotosynthesis, leaf::LeafPhotosynthesisParams, met::meteo, canRad)
  ISun   = CartesianIndices(canRad.absPAR_sunCab)
  IShade = CartesianIndices(canRad.absPAR_shadeCab)
  for i in IShade
    leaf.APAR = 1e6*canRad.absPAR_shadeCab[i]
    leaf.T    = canRad.T_shade[i]
    LeafPhotosynthesis!(mo, leaf, met)
    canRad.GPP_shade[i] = leaf.An;
    canRad.ϕ_shade[i]   = leaf.ϕs;
    canRad.gs_shade[i]  = leaf.gs;
    canRad.H_shade[i]   = leaf.H;
    canRad.LE_shade[i]  = leaf.LE;
    canRad.Cc_shade[i]  = leaf.Cc;
    canRad.NPQ_shade[i] = leaf.NPQ;
  end
  for i in ISun
    leaf.APAR = 1e6*canRad.absPAR_sunCab[i]
    leaf.T    = canRad.T_sun[i[3]]
    LeafPhotosynthesis!(mo, leaf, met)
    canRad.GPP_sun[i] = leaf.An;
    canRad.ϕ_sun[i]   = leaf.ϕs;
    canRad.gs_sun[i]  = leaf.gs;
    canRad.H_sun[i]   = leaf.H;
    canRad.LE_sun[i]  = leaf.LE;
    canRad.Cc_sun[i]  = leaf.Cc;
    canRad.NPQ_sun[i] = leaf.NPQ;
  end
end

=#