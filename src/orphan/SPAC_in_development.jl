###############################################################################
#
# Leaf boundary layer parameter set
#
###############################################################################
#= AbstractLeafBLParaSet type tree
---> LeafBLParaSetFixed
---> LeafBLParaSetGentine
=#

#=
"""
    AbstractLeafBLParaSet

Hierarchy of the `AbstractLeafBLParaSet`:
- [`LeafBLParaSetFixed`](@ref)
- [`LeafBLParaSetGentine`](@ref)
"""
abstract type AbstractLeafBLParaSet end




"""
    struct LeafBLParaSetFixed{FT<:AbstractFloat}

Leaf boundary layer with fixed resistance

# Fields
$(DocStringExtensions.FIELDS)
"""
struct LeafBLParaSetFixed{FT<:AbstractFloat} <: AbstractLeafBLParaSet
    "Resistance of the leaf boundary layer `[m² s mol⁻¹]`"
    ra::FT
end




"""
    struct LeafBLParaSetGentine

Gentine's boundary layer scheme that computes Monin Obhukow length and
    resistances across the leaf to Ca scale
"""
struct LeafBLParaSetGentine <: AbstractLeafBLParaSet end
=#





###############################################################################
#
# Grid information struct
# This struct should go to types.jl
#
###############################################################################
#=
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
    l.Rn      =   (1-l.α)*met.S_down +  2 * (met.L_down - l.ε* K_STEFAN *Tleaf^4); # 2 is for two sides of the leaves
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


###############################################################################
#
# Update leaf boundary layer conductance
#
###############################################################################



#=
function boundary_layer_resistance!(bl::LeafBLParaSetGentine, leaf::LeafParams,  envir::MeteoParams)
    # based on Monin-Obukhov Similiarity theory -> to be changed for LES
    # compute Obukhov length
    # iterate a few times

    # TODO ideally should not have any information about the leaves as it is the turbuence above canopy in log profile
    # first update roughness (if phenology is changing)
    # need to make this abstract as well.
    setRoughness!(leaf)

    if(leaf.height>envir.zscreen)
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
    Tv       = envir.t_air*(1.0+0.61*envir.e_air);
    ra_leaf  = 1/leaf.Cd / sqrt(envir.U/leaf.dleaf);
    #println("ra_leaf=",ra_leaf)
    ΔT   =   leaf.T - envir.t_air;
    VPD      =   max(leaf.esat-envir.e_air,1.0)
    ρd       =   envir.P_air/(R_D * envir.t_air)

    # TODO use ThermoDynamics to compute latent heat flux
    lv       =   latent_heat_vapor(leaf.T);
    L        =   envir.L; # initial Obukhov length
    counter = 1;
    #@show abs(1.0-Lold/L)
    while (counter<20 && abs(1.0-Lold/L)>1e-4) # 1% error
      #println("L=",L," ,Lold=",Lold)
      ra_m     =   max(1.0/(VON_KARMAN_CONST^2*envir.U) * ( log((envir.zscreen - leaf.d)/leaf.z0m) - ψ_m((envir.zscreen - leaf.d)/L,envir.stab_type_stable) + ψ_m(leaf.z0m/L,envir.stab_type_stable) ) * ( log((envir.zscreen - leaf.d)/leaf.z0m) - ψ_m((envir.zscreen - leaf.d)/L,envir.stab_type_stable) + ψ_h(leaf.z0h/L,envir.stab_type_stable) ), rmin) ;# momentum aerodynamic resistance
      ra_w     =   max(1.0/(VON_KARMAN_CONST^2*envir.U) * ( log((envir.zscreen - leaf.d)/leaf.z0m) - ψ_m((envir.zscreen - leaf.d)/L,envir.stab_type_stable) + ψ_m(leaf.z0m/L,envir.stab_type_stable) ) * ( log((envir.zscreen - leaf.d)/leaf.z0h) - ψ_h((envir.zscreen - leaf.d)/L,envir.stab_type_stable) + ψ_h(leaf.z0h/L,envir.stab_type_stable) ), rmin) ;# water aerodynamic resistance
      #println("ra_m=",ra_m)
      #@show ra_leaf
      #@show ra_w
      ram_full =   ra_leaf + ra_m;
      raw_full =   ra_leaf + ra_w;
      H        =   ρd * CP_D * ΔT / raw_full;
      rs_s_m   =   envir.g_m_s_to_mol_m2_s/leaf.gs;
      LE       =   WATER_AIR_MRATIO / envir.P_air * ρd * lv * VPD / (rs_s_m + raw_full)
      #@show LE
      ustar    =   sqrt(envir.U/ram_full);
      Hv_s     =   H + 0.61 * CP_D / lv * envir.t_air * LE
      Lold     =   L;
      L        =   - ustar^3*Tv/(GRAVITY * VON_KARMAN_CONST * Hv_s); # update Obukhov length
      #println("L=",L, " ra=",ra_w," (s/m), H=", H, " (W/m2), counter=", counter)
      counter = counter+1
    end

    # TODO Is the envir more a "global" variable set for the tree? Or a "local" variable for each layer/leaf?
    # If for the whole plant or layer, then envir should be updated only once?
    # save these values in leaf and flux structures
    envir.L   = L
    leaf.H  = H
    leaf.LE = LE
    envir.ustar = ustar
    envir.ra = raw_full
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
    leaf.z0m = 0.1*leaf.height;    # tree roughness (m)
    leaf.z0h = 0.1*leaf.height;    # tree roughness (m) - TODO should be changed later
    leaf.d   = 2/3*leaf.height;    # tree displacement height (m)
end
=#
