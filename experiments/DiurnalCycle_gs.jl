
#using Revise
#using PyPlot

using Land
using Land.PhysCon
using Land.WaterVapor
using Land.Leaf
#using IncGammaBeta
output_dir = joinpath(@__DIR__,"experiments/output")
mkpath(output_dir)

using Plots
using BenchmarkTools
using DifferentialEquations



# Create a leaf structure
l = leaf_params{Float32}();
# Create a Flux structure
f = fluxes{Float32}();
# Create a meteo structure
met = meteo{Float32}();

#met = meteo();

# initialize some reasonable values
f.Je   = 100; f.gbc  = 100; f.gbv  = 100; f.ceair= 1500; f.eair = 1500; f.APAR = 500; f.H=0;f.LE=0; # leaf should not have eair
l.Kn = 2.44; l.α=0.2; l.ε=0.98; l.LMA=100e-3; l.RWC=80/100;l.psi_l=-1e6;l.psi_l50 = -1e6;l.ck=3;met.zscreen = 2.0;
l.height   = 1.0; met.zscreen  = 2.0;
met.stab_type_stable = 2;
l.gstyp = 1;
l.dynamic_state=true


# A diurnal cycle for radiation and Tair
Deltat  = 60;
Samp    = 500; # W/m2 amplitude
Tmean   = 273.15+15;
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



plot(t/3600,Sdown_t,label="Sdown (W/m2)")
plot!(t/3600,10*(Tair_t-273.15*ones(size(Tair_t))),label="Tair (C)")

psi_s      = -1e6 ; # soil water potential (Pa)
U          =  1.0;
RH         =  65/100;
eps_air    =  0.75;
tspan      =  (0.0,Deltat);
N          =  length(Sdown_t);
mutable struct parameters_ode
    l::leaf_params;
    met::meteo;
    f::fluxes;
    psi_s;
end



# small time stepping
dt     = 0.1*60; # in s
T_t    = zeros(size(Sdown_t));
psil_t = zeros(size(Sdown_t));
Cc_t   = zeros(size(Sdown_t));
gs_t = zeros(size(Sdown_t));
Rn_t   = zeros(size(Sdown_t));
GPP_t  = zeros(size(Sdown_t));
GPP_diffusion_t  = zeros(size(Sdown_t));
LUE_t  = zeros(size(Sdown_t));
H_t    = zeros(size(Sdown_t));
LE_t   = zeros(size(Sdown_t));
rs_t   = zeros(size(Sdown_t));
ra_t   = zeros(size(Sdown_t));

function f_ode!(du,u,p,t) # p are parameters
    du .= LeafEnergyWaterBalance(u[1], u[2], u[3], p.met, p.l, p.f, p.psi_s);
    #println("du_inside = $(du), u_inside = $(u)")
end

Sdown_t[300:400].=600;
Sdown_t[800:1000].=1400;
let
    # initial conditions
    met.T_air  = Tair_t[1];
    l.T        = met.T_air;
    l.psi_l    = psi_s;
    l.Cc       = 0.6*met.Ca;
    #println("Tair1=",met.T_air," Tleaf1=",l.T," psi_leaf1=",l.psi_l)

    for i=1:N
        met.S_down = Sdown_t[i];
        met.L_down = eps_air*physcon.σ*(Tair_t[i])^4;
        met.T_air  = Tair_t[i];
        met.e_air  = RH*SatVap(Tair_t[i])[1];
        met.PAR    = 45/100*physcon.Wtoμmole_s*Sdown_t[i];
        f.APAR     = met.PAR;
        met.U      = U;
        met.Ca     = 400.0; #ppm
        f.APAR     = met.PAR;
        #println("Tair=",met.T_air," Tleaf=",l.T," psi_leaf=",l.psi_l)
        for j=1:trunc(Deltat/dt)
            u    = [l.T;l.psi_l;l.Cc;l.gs];
            p    = parameters_ode(l,met,f,psi_s);
            #(p.met, p.l,  p.psi_s, p.U)    = [l;met;psi_s;U];
            #prob = ODEProblem(f_ode!,u0,tspan,p);
            #du   = zeros(size(u));
            #f_ode!(du,u,p,t,dummy);
            #println("du_outside = $(du), u_outside = $(u)")
            # Rn_t[i] = dummy[3]; H_t[i] = dummy[4]; LE_t[i] = dummy[5];


            du   = zeros(size(u));
            f_ode!(du,u,p,t);
            # just remove Cc part here.
            du[3]=0
            #@show l.gs
            (l.T,l.psi_l,test, l.gs) = du*dt+u;
            #@show l.gs
#            u0   = [l.T;l.psi_l;l.Cc];
#            prob = ODEProblem(f_ode!,u0,tspan,p);
#            sol  = solve(prob);
#            # save values
#            met = p.met;
#            f   = p.f;
#            l   = p.l;
#            (l.T,l.psi_l,l.Cc) = sol[1:3,end];

            #println("Cc_out=",l.Cc)
            Rn_t[i] = p.f.Rn; H_t[i] = p.f.H; LE_t[i] = p.f.LE;
            rs_t[i] = 1.0/(p.l.gs/p.f.g_m_s_to_mol_m2_s);
            ra_t[i] = p.f.ra;
            GPP_t[i]= p.f.An_biochemistry;
            GPP_diffusion_t[i]= p.f.An_diffusion;
            LUE_t[i]= p.f.An_biochemistry/f.APAR;

            #T_old   = l.T;
            #(l.T,l.psi_l) = du*dt+u;
            T_t[i]    = l.T ; #  = T_old;
            psil_t[i] = l.psi_l;
            Cc_t[i]   = l.Cc;
            gs_t[i]   = l.gs
            if(abs(H_t[i])>500)
                println("index ($i) ($j)")
            end


#             (dumb,l.psi_l) = du*dt+u;
#             T_t[i]  = dumb ; #  = T_old;


        end
    end




end


# met.L = 1e6;
# setra!(l, f, met) ;
# log((met.zscreen - l.d)/l.z0m)
# l.height



plot(t/3600, T_t-273.15*ones(size(T_t)),xlabel = "t (hr)",ylabel = "T (C)",label="T (C)",ylim=0:100)

plot(t/3600,psil_t/1e6,xlabel = "t (hr)",ylabel = "|psi_l (MPa)|",label="psi_l",ylim=0:100)

savefig(joinpath(output_dir, "T_psi_diurnal_t_gs.png"))

plot(t/3600,  Rn_t,xlabel = "t (hr)",ylabel = "Rn (W/m2)",label="Rn")
plot!(t/3600, H_t,xlabel  = "t (hr)",ylabel = "H (W/m2)" ,label="H")
plot!(t/3600, LE_t,xlabel = "t (hr)",ylabel = "LE (W/m2)",label="LE")

savefig(joinpath(output_dir, "Fluxes_diurnal_t_gs.png"))

plot(t/3600,  rs_t,xlabel = "t (hr)",ylabel = "rs (s/m)",label="rs")
plot!(t/3600, ra_t,xlabel  = "t (hr)",ylabel = "ra (s/m)" ,label="ra")

savefig(joinpath(output_dir, "resistances_diurnal_t_gs.png"))

plot(t/3600, GPP_t,xlabel = "t (hr)",ylabel = "GPP (micomles/m2/s)",label="GPP",ylim=0:100)

savefig(joinpath(output_dir, "GPP_diurnal_t_gs.png"))

plot(t/3600, LUE_t,xlabel = "t (hr)",ylabel = "LUE (-)",label="GPP",ylim=0:100)

savefig(joinpath(output_dir, "LUE_diurnal_t_gs.png"))

plot(t/3600, gs_t,xlabel = "t (hr)",ylabel = "gs (-)",label="gs",ylim=0:100)

savefig(joinpath(output_dir, "Gs_diurnal_gs.png"))

plot(t/3600, Cc_t,xlabel = "t (hr)",ylabel = "Cc (-)",label="Cc",ylim=0:100)

savefig(joinpath(output_dir, "Cc_diurnal_gs.png"))

plot(t/3600, T_t,xlabel = "t (hr)",ylabel = "T_leaf (-)",label="Tleaf",ylim=0:100)

savefig(joinpath(output_dir, "Tleaf_diurnal_gs.png"))


