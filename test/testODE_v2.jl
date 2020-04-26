using DifferentialEquations
using Interpolations
using Plots

@testset "Test ODE, version 2" begin

   theta_res = 0.075;     # Residual water content
   theta_sat = 0.287;     # Volumetric water content at saturation
   vg_alpha = 0.027;      # Inverse of the air entry potential [/cm]
   vg_n = 3.96;           # Pore-size distribution index
   vg_m = 1;              # Exponent
   Ksat = 34 / 3600;      # Hydraulic conductivity at saturation [cm/s]

   params = [theta_res theta_sat vg_alpha vg_n vg_m Ksat];
   # Van Genuchten Function

   function vG(ψ; params=params)
   # ----------------------------------
   # van Genuchten [1980] relationships
   # ----------------------------------
   # --- Soil parameters
   θ_res = params[1];   # Residual water content
   θ_sat = params[2];   # Volumetric water content at saturation
   α = params[3];       # Inverse of the air entry potential
   n = params[4];           # Pore-size distribution index
   m = params[5];           # Exponent
   Ksat = params[6];        # Hydraulic conductivity at saturation

   # --- Effective saturation [Se] for specified matric potential [psi]
   if ψ < 0
      Se = (1 + (α * abs(ψ))^n)^-m
   else
      Se = 1
   end

   # --- Volumetric soil moisture [theta] for specified matric potential [psi]
   θ = θ_res + (θ_sat - θ_res) * Se

   # --- Hydraulic conductivity [K] for specified matric potential [psi]
   if Se < 1
      #K = Ksat * sqrt(Se) * (1 - (1 - Se^(1/m))^m)^2
      K = Ksat * 1.175e6 / (1.175e6 + abs(ψ)^4.74);
   else
      K = Ksat
   end

   # --- Specific moisture capacity [cap] for specified matric potential [psi]
   if (ψ < 0)
      num = α * m * n * (θ_sat - θ_res) * (α * abs(ψ))^(n-1)
      den =  (1 + (α * abs(ψ))^n)^(m+1)
      cap = num / den
   else
      cap = 0.0
   end

   return θ, K, cap
   end

   # Discretization (should be flexible in delta-z, need to try out later)
   N = 150
   # Let's do layers
   z_bnd = zeros(N+1);
   for i=2:length(z_bnd)
      z_bnd[i] = z_bnd[i-1]-i-1;
   end
   # Layer thickness:
   Δz = z_bnd[1:end-1]-z_bnd[2:end];
   # Layer center:
   z = 0.5*(z_bnd[1:end-1]+z_bnd[2:end]);
   # Half Layer plus/minus thickness
   Δz_half = similar(z_bnd);
   Δz_half[1]=z_bnd[1]-z[1];
   for i = 2:N
      Δz_half[i] = z[i-1]-z[i];
   end
   Δz_half[N+1] = z[N]-z_bnd[N+1];

   #Δx[1:10]=-1
   # Set prior Ψ
   Ψ_0 = zeros(N);
   θ_0 = zeros(N);
   θ_0[:] .= 0.1;

   Ψ_0[:] .= -100;
   #Ψ_0[50:end]=-
   #Ψ_0[50:end] .= -

   function matric_potential(θ::Number; params=params)
   # --- Calculate psi for a given theta
      θ_res = params[1];    # Residual water content
      θ_sat = params[2];    # Volumetric water content at saturation
      α = params[3];        # Inverse of the air entry potential
      n = params[4];            # Pore-size distribution index
      m = params[5];            # Exponent
      if θ<=θ_res
         θ=max(θ_res,θ)+1e-10
      elseif θ>=θ_sat
         θ=min(θ_sat,θ)-1e-10
      end

      Se = (θ  - θ_res) / (θ_sat - θ_res);
      psi = -((Se^(-1/m) - 1)^(1/n)) / α;
   end




   # Still too many allocations here but we can fine-tune later
   #This basically solves Eq 8.24 in Gordon's book
   function soil_water(du,u,p,t)
      re = vG.(u, params=p)
      # need interpolation to get to half levels
      #theta = [x[1] for x in re]
      KK = [x[2] for x in re]
      C = [x[3] for x in re]
      #KK[end]=0.0
      for i=1:length(u)
         # Top layer
         if i==1
            Kp = 0.5*(KK[i]+KK[i+1])
            Km = 0.0037
            # Just as test, boundary condition of Psi=-20 here.
            du[i] = (Km/Δz_half[i]*(-20-u[i])-Kp/Δz_half[i+1]*(u[i]-u[i+1])+Km-Kp)/Δz[i]/C[i]
            # Just as test, set flux boundary here.


         # Bottom Layer
         elseif i==length(u)
            Km = 0.5*(KK[i]+KK[i-1])
            # Free drainage at the bottom (prefer to model groundwater later):
            #du[i] = (Km/Δz[i]*(u[i-1]-u[i])+Km-KK[i])/Δz[i]/C[i]
            du[i] = (Km/Δz[i]*(u[i-1]-u[i])+Km-KK[i])/Δz[i]/C[i]

         else
            Km = 0.5*(KK[i]+KK[i-1])
            Kp = 0.5*(KK[i]+KK[i+1])
            du[i] = (Km/Δz_half[i]*(u[i-1]-u[i])-Kp/Δz_half[i+1]*(u[i]-u[i+1])+Km-Kp)/Δz[i]/C[i]
         end
      end
   end


   tspan = (0.0,20*3600.0)
   p = params
   prob = ODEProblem(soil_water,Ψ_0,tspan,p)
   alg = ImplicitEuler()
   # Save every 10min
   @time sol = solve(prob, alg,saveat=10*60, dtmax=5);

   plot(sol[:,1:end], z/100)

end
