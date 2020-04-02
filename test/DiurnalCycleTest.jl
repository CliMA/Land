
# Add PATH
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

#using PyPlot
using Plots
using BenchmarkTools
using DifferentialEquations

using LSM.PhysCon
using LSM.WaterVaporMod
using LSM.LeafPhotosynthesisMod

output_dir = joinpath(@__DIR__,"..","output")
mkpath(output_dir)

# function lorenz!(du,u,p,t)
#  du[1] = 10.0*(u[2]-u[1]);
#  du[2] = u[1]*(28.0-u[3]) - u[2];
#  du[3] = u[1]*u[2] - (8/3)*u[3];
# end
# u0 = [1.0;0.0;0.0];
# tspan = (0.0,100.0);
# prob = ODEProblem(lorenz!,u0,tspan);
# sol = solve(prob);

# plot(sol,vars=(1,2,3))

# Create a leaf structure
l = leaf_params{Float32}();
# Create a Flux structure
f = LeafPhotosynthesisMod.fluxes{Float32}();
# Create a meteo structure
met = meteo{Float32}();
#met = meteo();

# initialize some reasonable values
f.je   = 100; f.gbc  = 100; f.gbv  = 100; f.ceair= 1500; f.eair = 1500; f.APAR = 500; # leaf should not have eair
l.Kn = 2.44; l.α=0.2; l.ε=0.98; l.LMA=40e-3; l.c_leaf=50/100*4184;l.ra=50;l.psi_l=-1e6;l.psi_l50 = -2e6;l.ck=3;

# A diurnal cycle for radiation and Tair
Deltat  = 60;
Samp    = 500; # W/m2 amplitude
Tmean   = 273.15+22;
DeltaT  = 3;
omega   = 2*π/(24*3600);
t       = range(0, stop=24*3600, step=Deltat); # diurnal cycle in seconds
#print(t)
phi_t   = omega*t-π*ones(size(t))/2;
Sdown_t = zeros(size(t));
Tair_t  = zeros(size(t));
zeros_t = zeros(size(t));
for i = 1:length(Sdown_t)
    Sdown_t[i] = Samp*max( sin(phi_t[i]),zeros_t[i] );#max(sin(phase[i]),zeros(size(t[i])));
    Tair_t[i]  = Tmean + DeltaT*sin(phi_t[i]-π/3);
end
#print(typeof(phase[1]));
#print(typeof(Sdown_t[1]));
#print(Sdown_t);
#clf();
#fig = figure("Diurnal cycle", figsize=(10,5));
plot(t/3600,Sdown_t)
plot!(t/3600,10*(Tair_t-273.15*ones(size(Tair_t))))
#    title = "Shortwave Incoming radiation (W/m^2)");#,
#    xlabel = "Hours",
#    ylabel = "S_{down} (W/m^2)")
#plot(t/3600,Tair_t-273.15*ones(size(Tair_t)),
#    title = "Temperature (C)",
#    xlabel = "Hours",
#    ylabel = "Temperature (C)")

psi_s      = -0.3e6 ; # soil water potential (Pa)
U          =  1;
RH         =  65/100;
eps_air    =  0.75;
tspan      =  (0.0,Deltat);
N          =  length(Sdown_t);
mutable struct parameters_ode
    l::leaf_params;
    met::meteo;
    psi_s;
    U;
end

# small time stepping
dt = 0.1*60; # in s
T_t    = zeros(size(Sdown_t));
psil_t = zeros(size(Sdown_t));


function f_ode!(du,u,p::parameters_ode,t) # p are parameters
    du = LeafEnergyWaterBalance(u[1], u[2], p.met, p.l,  p.psi_s, p.U);
    println("du_inside = $(du), u_inside = $(u)")
end

let

    for i=1:1;#N
        met.S_down = Sdown_t[i];
        met.L_down = eps_air*physcon.σ*(Tair_t[i])^4;
        met.T_air  = Tair_t[i];
        met.ea_air = RH*SatVap(Tair_t[i])[1];
        met.PAR    = 45/100*physcon.Wtoμmole_s*Sdown_t[i];
        met.Cs     = 400.0; #ppm
        for j=1:trunc(Deltat/dt)
            u    = [l.T;l.psi_l];
            p    = parameters_ode(l,met,psi_s,U);
            #(p.met, p.l,  p.psi_s, p.U)    = [l;met;psi_s;U];
            #prob = ODEProblem(f_ode!,u0,tspan,p);
            du   = zeros(size(u));
            f_ode!(du,u,p,t)
            println("du_outside = $(du), u_outside = $(u)")
            (l.T,l.psi_l) = du*dt+u;
            T_t[i] = l.T; psil_t[i]=l.psi_l;
        end
    end

    plot(T_t)
    savefig(joinpath(output_dir, "T_t.png"))


    l.T        = Tair_t[1]; # initialize temperature of the leaf
    l.psi_l    = psi_s;

    for i=1:1;N
        met.S_down = Sdown_t[i];
        met.L_down = eps_air*physcon.σ*(Tair_t[i])^4;
        met.T_air  = Tair_t[i];
        met.ea_air = RH*SatVap(Tair_t[i])[1];
        u0   = [l.T;l.psi_l];
        p    = parameters_ode(l,met,psi_s,U);
        #(p.met, p.l,  p.psi_s, p.U)    = [l;met;psi_s;U];
        prob = ODEProblem(f_ode!,u0,tspan,p);
        sol  = solve(prob);
        (l.T,l.psi_l) = sol[:,end];
        T_t[i] = l.T; psil_t[i]=l.psi_l;
        #print((l.T,l.psi_l) )
        #dT_dt,dH2Ol_dt = LeafEnergyWaterBalance(met, l, psi_s);
        #l.T        = l.T + Deltat*dT_dt;
    end

    plot(T_t)
    savefig(joinpath(output_dir, "T_t_final.png"))

end

# LeafEnergyWaterBalance(met, l, psi_s)
