
#using Revise

# Add PATH
#push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

#using PyPlot
using Plots
using BenchmarkTools
using DifferentialEquations

using Land
using Land.PhysCon
using Land.WaterVapor
using Land.Leaf

output_dir = joinpath(@__DIR__,"output")
#mkpath(output_dir)

using DelimitedFiles
PAR  = readdlm("/Users/cfranken/PAR.dat");
size(PAR)
# use just a subset 
#PAR = PAR[1:1000,:];

tmax = length(PAR) # Time in seconds here
t = range(1,tmax,step=1)
l = leaf_params{Float32}()
f = fluxes{Float32}()
met = meteo{Float32}()

mutable struct parameters_ode
    l::leaf_params;
    met::meteo;
    f::fluxes;
    psi_s;
end

function f_ode!(du,u,p,t) # p are parameters
    du .= LeafEnergyWaterBalance(u[1], u[2], u[3], p.met, p.l, p.f, p.psi_s)[1:3];
    #println("du_inside = $(du), u_inside = $(u)")
end


f.Je   = 100; f.gbc  = 100; f.gbv  = 100; f.ceair= 1500; f.eair = 1500;  f.H=0;f.LE=0; # leaf should not have eair
l.Kn = 2.44; l.α=0.2; l.ε=0.98; l.LMA=100e-3; l.RWC=80/100;l.psi_l=-1e6;l.psi_l50 = -2e6;l.ck=3;met.zscreen = 2.0;
l.height   = 1.0; 
met.stab_type_stable = 2;
psi_s      = -0.3e6 ; # soil water potential (Pa)
U          =  1.0;
RH         =  65/100;
eps_air    =  0.75;

l.dynamic_state = true
l.gstyp = 1 # Ball Berry
l.g1_BB=4
met.Ca = 400
l.Vcmax25 = 90
l.Jmax25 = 90*1.9
uu = zeros(tmax,16)
uu[1,2] = 0.07
uu[1,1] = 0.5
f.APAR = 50
l.Rdleaf = 1.0
met.T_air = 298;
l.T = met.T_air + 0.5;
l.psi_l = psi_s;
l.Cc = 0.6*400;
#apar[1:500]=0
#apar[501:100]=0
tspan      =  (0.0,1.0);


dt     = 0.1*60
for c = 2:1:tmax
    #println(apar[c])
       met.S_down = PAR[c]/(45/100*physcon.Wtoμmole_s);
       #println(met.S_down)
       met.T_air  = 298.0;
       met.L_down = eps_air*physcon.σ*(met.T_air)^4;
       met.e_air  = RH*SatVap(met.T_air)[1];
       met.PAR    = PAR[c];
       met.U      = 1.0;
    
       f.APAR = PAR[c]*0.5 
       f.APAR = max(10,f.APAR)
       #l.gs = uu[c-1,2]
       l.Kn = uu[c-1,1]
       uu[c-1,3]=l.ϕs
       uu[c-1,4]=l.Cc
       uu[c-1,5]=f.φ
       uu[c-1,6]=l.Kp
       uu[c-1,7]=f.APAR
       uu[c-1,8]=f.An_biochemistry
       uu[c-1,9]=l.Fm′
       uu[c-1,10]=l.T
       uu[c-1,11]=f.Rn
       uu[c-1,12]=f.H
       uu[c-1,13]=f.LE
       uu[c-1,14]=l.Cc
       uu[c-1,15]=f.An_diffusion
       uu[c-1,16]=l.gs
#        println("gs out=",l.gs)
       #l.Kp = 4
       u    = [l.T;l.psi_l;l.Cc];
       p    = parameters_ode(l,met,f,psi_s);
       du   = zeros(size(u));
       f_ode!(du,u,p,t);
       # just remove Cc part here.
       #du[3]=0
            #@show l.gs
        (l.T,l.psi_l,l.Cc) = du*dt+u;
        #(p.met, p.l,  p.psi_s, p.U)    = [l;met;psi_s;U];
        #prob = ODEProblem(f_ode!,u0,tspan,p);
        #du   = zeros(size(u));
        #f_ode!(du,u,p,t,dummy);
        #println("du_outside = $(du), u_outside = $(u)")
        # Rn_t[i] = dummy[3]; H_t[i] = dummy[4]; LE_t[i] = dummy[5];


#             du   = zeros(size(u));
#             f_ode!(du,u,p,t);
#             (l.T,l.psi_l,l.Cc) = du*dt+u;

       
            
    
       #LeafPhotosynthesis!(f,l,met)
    if (l.Kn_ss-l.Kn) > 0
        tau_k = 1
    else
        tau_k = 1
    end
       uu[c,1] = uu[c-1,1]+(l.Kn_ss-l.Kn)/tau_k*1/60
       uu[c,2] = uu[c-1,2]+(l.gs_ss-l.gs)/20*1/60
    #println(l.Kn_ss, " ",  l.Ci, " ", f.φ, " ", l.gs_ss, " ", f.APAR)
end

plot(t/3600, uu[:,8],label="GPP")
plot!(t/3600, PAR/10,label="PAR/10")

plot(t/3600, uu[:,16],label="g_leaf")

plot(t/3600, uu[:,8],label="GPP biochemistry")
plot!(t/3600, uu[:,15],label="GPP diffusion")

plot(t/3600, uu[:,11],label="Rn",ylim=-50:500)
plot!(t/3600, uu[:,12],label="H",ylim=-50:500)
plot!(t/3600, uu[:,13],label="LE",ylim=-50:500)

plot(t/3600, uu[:,14]/met.Ca,label="Cc/Ca")










