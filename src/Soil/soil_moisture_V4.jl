module soil_moisture_V4
using Parameters

export soil_struct,
       van_Genuchten,
       Campbell,
       tridiagonal_solver,
       matric_potential,
       root_uptake,
       predictor_corrector,
       compute_grid_settings

# --- Define soil layers

# Number of soil layers
@with_kw mutable struct soil_struct{}
   nsoi            # Number of soil layers
   dz              = ones(nsoi)       # Soil layer thickness (cm)
   z_plus_onehalf  = zeros(nsoi) # Soil depth [cm] at i+1/2 interface between layers i & i+1 [negative distance from surface]
   z               = zeros(nsoi) # Soil depth [cm] at center of layer i [negative distance from surface]
   dz_plus_onehalf = zeros(nsoi) # Thickness between between z[i] & z[i+1]
   functions       = "aa"
   theta           = zeros(nsoi) # Soil Moisture
   psi             = zeros(nsoi) # Matric potential
   theta0          = 0.0  # Initial SM
   psi0            = 0.0  # initial matric potential
   K               = zeros(nsoi) # Hydraulic conductivity (cm H2O/s)
   cap             = zeros(nsoi) # Specific moisture capacity (/cm)
   Q0              = 0.0 # Infiltration flux (cm H2O/s)
   QN              = 0.0 # Drainage flux (cm H2O/s)
   dtheta          = 0.0 # Change in soil moisture (cm H2O)
   err             = 0.0 # Water balance error (cm H2O)

   # Some root distribution functions
   ssflag          = 1 # perform root uptake else do not
   bi              = 0.98
   fz              = zeros(nsoi)  # root distribution function
   psidry          = -20
   psiopt          = -60
   sink            = zeros(nsoi)
   beta            = zeros(nsoi)  # soil wetness factor

end


# Van Genuchten Function
function van_Genuchten(params, psi)

# ----------------------------------
# van Genuchten [1980] relationships
# ----------------------------------

# --- Soil parameters

theta_res = params[1];   # Residual water content
theta_sat = params[2];   # Volumetric water content at saturation
alpha = params[3];       # Inverse of the air entry potential
n = params[4];           # Pore-size distribution index
m = params[5];           # Exponent
Ksat = params[6];        # Hydraulic conductivity at saturation
ityp = params[7];        # Soil texture flag

# --- Effective saturation [Se] for specified matric potential [psi]

if (psi <= 0)
   Se = (1 + (alpha * abs(psi))^n)^-m
else()
   Se = 1
end

# --- Volumetric soil moisture [theta] for specified matric potential [psi]

theta = theta_res + (theta_sat - theta_res) * Se

# --- Hydraulic conductivity [K] for specified matric potential [psi]

if (Se <= 1)
   K = Ksat * sqrt(Se) * (1 - (1 - Se^(1/m))^m)^2

   # Special case for Haverkamp et al. (1977) sand [ityp = 1] & Yolo light clay [ityp = 2]

   if (ityp == 1)
      K = Ksat * 1.175e6 / (1.175e6 + abs(psi)^4.74)
   end
   if (ityp == 2)
      K = Ksat * 124.6/ (124.6 + abs(psi)^1.77)
   end

else()

   K = Ksat

end

# --- Specific moisture capacity [cap] for specified matric potential [psi]

if (psi <= 0)
   num = alpha * m * n * (theta_sat - theta_res) * (alpha * abs(psi))^(n-1)
   den =  (1 + (alpha * abs(psi))^n)^(m+1)
   cap = num / den
else()
   cap = 0
end

return theta, K, cap

end # function

# Van Genuchten Function
function vG(ψ, params=params)
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
if (ψ <= 0)
   Se = (1 + (α * abs(ψ))^n)^-m
else()
   Se = 1
end

# --- Volumetric soil moisture [theta] for specified matric potential [psi]
θ = θ_res + (θ_sat - θ_res) * Se

# --- Hydraulic conductivity [K] for specified matric potential [psi]
if Se <= 1
   K = Ksat * sqrt(Se) * (1 - (1 - Se^(1/m))^m)^2
else
   K = Ksat
end

# --- Specific moisture capacity [cap] for specified matric potential [psi]
if (ψ <= 0)
   num = α * m * n * (θ_sat - θ_res) * (α * abs(ψ))^(n-1)
   den =  (1 + (α * abs(ψ))^n)^(m+1)
   cap = num / den
else
   cap = 0
end

return θ, K, cap
end # function

function Campbell(params, psi)

# -----------------------------
# Campbell (1974) relationships
# -----------------------------

# --- Soil parameters

theta_sat = params[1];    # Volumetric water content at saturation
psi_sat = params[2];      # Matric potential at saturation
b = params[3];            # Exponent
Ksat = params[4];         # Hydraulic conductivity at saturation

# --- Volumetric soil moisture [theta] for specified matric potential [psi]

if (psi <= psi_sat)
   theta = theta_sat * (psi / psi_sat)^(-1/b)
else()
   theta = theta_sat
end

# --- Hydraulic conductivity [K] for specified matric potential [psi]

if (psi <= psi_sat)
   K = Ksat * (theta / theta_sat)^(2*b+3)
else()
   K = Ksat
end

# --- Specific moisture capacity [cap] for specified matric potential [psi]

if (psi <= psi_sat)
   cap = -theta_sat / (b * psi_sat) * (psi / psi_sat)^(-1/b-1)
else()
   cap = 0
end

return theta, K, cap

end # Function

# Other Functions
function tridiagonal_solver(a, b, c, d, n)

# Solve for U given the set of equations R * U = D; where U is a vector
# of length N; D is a vector of length N; & R is an N x N tridiagonal
# matrix defined by the vectors A, B, C each of length N. A[1] &
# C[N] are undefined and are not referenced.
#
#     |B[1] C[1] ...  ...  ...                     |
#     |A[2] B[2] C[2] ...  ...                     |
# R = |     A[3] B[3] C[3] ...                     |
#     |                    ... A[N-1] B[N-1] C[N-1]|
#     |                    ... ...    A[N]   B[N]  |
#
# The system of equations is written as:
#
#    A_i * U_i-1 + B_i * U_i + C_i * U_i+1 = D_i
#
# for i = 1 to N. The solution is found by rewriting the
# equations so that:
#
#    U_i = F_i - E_i * U_i+1

# --- Forward sweep [1 -> N] to get E & F
# Initialize E and F
e = copy(a)*0.0;
f = copy(a)*0.0;

e[1] = c[1] / b[1]


for i = 2: 1: n-1
   e[i] = c[i] / (b[i] - a[i] * e[i-1])
end


f[1] = d[1] / b[1]

for i = 2: 1: n
   f[i] = (d[i] - a[i] * f[i-1]) / (b[i] - a[i] * e[i-1])
end


# --- Backward substitution [N -> 1] to solve for U
u = zeros(n);
u[n] = f[n]

for i = n-1: -1: 1
   u[i] = f[i] - e[i] * u[i+1]
end

return u

end # function

function matric_potential(type, params, theta)

# --- Calculate psi for a given theta

if type == "van_Genuchten"

   theta_res = params[1];    # Residual water content
   theta_sat = params[2];    # Volumetric water content at saturation
   alpha = params[3];        # Inverse of the air entry potential
   n = params[4];            # Pore-size distribution index
   m = params[5];            # Exponent

   Se = (theta - theta_res) / (theta_sat - theta_res);
   psi = -((Se^(-1/m) - 1)^(1/n)) / alpha;

elseif type == "Campbell"

   theta_sat = params[1];    # Volumetric water content at saturation
   psi_sat = params[2];      # Matric potential at saturation
   b = params[3];            # Exponent

   psi = psi_sat * (theta / theta_sat)^-b;

end

return psi
end # function

function root_uptake(soil::soil_struct,ET)

   #Compute the soil sink terms

   densum = 0

   for i = 1:soil.nsoi
      if soil.psi[i] > soil.psidry
         soil.beta[i] = (soil.psi[i]-soil.psidry)/(soil.psiopt-soil.psidry)
      elseif soil.psi[i] <= soil.psidry
         soil.beta[i] = 0.0
      end

      densum = densum + soil.fz[i]*soil.beta[i]

   end

   # compute the sink terms

   for i = 1:soil.nsoi
      soil.sink[i] = ET*soil.fz[i]*soil.beta[i]/densum
   end

   return soil.sink

end # funnction


# Predictor Corrector Function
function predictor_corrector(soil::soil_struct, params, ET::Float64, dt)

# -------------------------------------------------------------
# Use predictor-corrector method to solve the Richards equation
# -------------------------------------------------------------

# Input
# dt                   ! Time step [s]
# soil.nsoi            ! Number of soil layers
# soil.functions       ! van Genuchten | Campbell relationships
# soil.dz_plus_onehalf ! Thickness between between z[i] & z[i+1] (cm)
# soil.dz              ! Soil layer thickness [cm]
# soil.psi0            ! Soil surface matric potential boundary condition [cm]
#
# Input/output
# soil.theta           ! Volumetric soil moisture
# soil.psi             ! Matric potential [cm]
#
# Output
# soil.K               ! Hydraulic conductivity [cm H2O/s]
# soil.cap             ! Specific moisture capacity [/cm]
# soil.Q0              ! Infiltration flux [cm H2O/s]
# soil.QN              ! Drainage flux [cm H2O/s]
# soil.dtheta          ! Change in soil moisture [cm H2O]
# soil.err             ! Water balance error (cm H2O)

# --- Save current soil moisture & matric potential for time n

theta_n = copy(soil.theta)
psi_n   = copy(soil.psi)

# --- Predictor step using implict solution for time n+1/2

# Hydraulic properties for current psi:
# theta - volumetric soil moisture
# K     - hydraulic conductivity
# cap   - specific moisture capacity

for i = 1:soil.nsoi
   if soil.functions == "van_Genuchten"
      soil.theta[i], soil.K[i], soil.cap[i] = van_Genuchten(params, soil.psi[i])
   elseif soil.functions == "Campbell"
      soil.theta[i], soil.K[i], soil.cap[i] = Campbell(params, soil.psi[i])
    end
end

# Hydraulic conductivity at i+1/2 interface between layers i & i+1 is the arithmetic mean()
K_plus_onehalf = soil.theta*0
for i = 1:soil.nsoi-1
   K_plus_onehalf[i] = 0.5 * (soil.K[i] + soil.K[i+1])
end

# Hydraulic conductivity at i=1/2 between surface (i=0) & first layer i=1

K_onehalf = soil.K[1]

# dz at i=1/2 between surface (i=0) & first layer i=1

dz_onehalf = 0.5 * soil.dz[1]

# Compute the sink terms
if soil.ssflag == 1

   soil.sink = root_uptake(soil,ET)
else

   soil.sink = soil.sink.*0
end



# Terms for tridiagonal matrix
a = zeros(soil.nsoi);
b = zeros(soil.nsoi);
c = zeros(soil.nsoi);
d = zeros(soil.nsoi);


i = 1
a[i] = 0
c[i] = -K_plus_onehalf[i] / soil.dz_plus_onehalf[i]
b[i] = soil.cap[i] * soil.dz[i] / (0.5 * dt) + K_onehalf / dz_onehalf - c[i]
d[i] = soil.cap[i] * soil.dz[i] / (0.5 * dt) * soil.psi[i] + K_onehalf / dz_onehalf * soil.psi0 + K_onehalf - K_plus_onehalf[i] - soil.sink[i]

for i = 2:soil.nsoi-1
   a[i] = -K_plus_onehalf[i-1] / soil.dz_plus_onehalf[i-1]
   c[i] = -K_plus_onehalf[i] / soil.dz_plus_onehalf[i]
   b[i] = soil.cap[i] * soil.dz[i] / (0.5 * dt) - a[i] - c[i]
   d[i] = soil.cap[i] * soil.dz[i] / (0.5 * dt) * soil.psi[i] + K_plus_onehalf[i-1] - K_plus_onehalf[i] - soil.sink[i]
end

i = soil.nsoi
a[i] = -K_plus_onehalf[i-1] / soil.dz_plus_onehalf[i-1]
c[i] = 0
b[i] = soil.cap[i] * soil.dz[i] / (0.5 * dt) - a[i] - c[i]
d[i] = soil.cap[i] * soil.dz[i] / (0.5 * dt) * soil.psi[i] + K_plus_onehalf[i-1] - soil.K[i] - soil.sink[i]

# Solve for psi at n+1/2

psi_pred = tridiagonal_solver(a, b, c, d, soil.nsoi)
# --- Corrector step using Crank-Nicolson solution for time n+1

# Hydraulic properties for psi_pred

for i = 1:soil.nsoi
   if soil.functions == "van_Genuchten"
      soil.theta[i], soil.K[i], soil.cap[i] = van_Genuchten(params, psi_pred[i])
   elseif  soil.functions == "Campbell"
      soil.theta[i], soil.K[i], soil.cap[i] = Campbell(params, psi_pred[i])
    end
end

# Hydraulic conductivity at i+1/2 interface between layers i & i+1

for i = 1:soil.nsoi-1
   K_plus_onehalf[i] = 0.5 * (soil.K[i] + soil.K[i+1])
end

# Hydraulic conductivity at i=1/2 between surface (i=0) & first layer i=1

K_onehalf = soil.K[1]

# dz at i=1/2 between surface (i=0) & first layer i=1

dz_onehalf = 0.5 * soil.dz[1]

# Compute the sink terms
if soil.ssflag == 1

   soil.sink = root_uptake(soil,ET)
else

   soil.sink = soil.sink.*0
end


# Terms for tridiagonal matrix

i = 1
a[i] = 0.0
c[i] = -K_plus_onehalf[i] / (2.0 * soil.dz_plus_onehalf[i])
b[i] = soil.cap[i] * soil.dz[i] / dt  + K_onehalf / (2.0 * dz_onehalf) - c[i]
d[i] = soil.cap[i] * soil.dz[i] / dt * soil.psi[i] + K_onehalf / (2.0 * dz_onehalf) * soil.psi0 + K_onehalf / (2.0 * dz_onehalf) * (soil.psi0 - soil.psi[i]) + c[i] * (soil.psi[i] - soil.psi[i+1]) + K_onehalf - K_plus_onehalf[i] - soil.sink[i]


for i = 2:soil.nsoi-1
   a[i] = -K_plus_onehalf[i-1] / (2.0 * soil.dz_plus_onehalf[i-1])
   c[i] = -K_plus_onehalf[i] / (2.0 * soil.dz_plus_onehalf[i])
   b[i] = soil.cap[i] * soil.dz[i] / dt - a[i] - c[i]
   d[i] = soil.cap[i] * soil.dz[i] / dt * soil.psi[i] - a[i] * (soil.psi[i-1] - soil.psi[i]) + c[i] * (soil.psi[i] - soil.psi[i+1]) + K_plus_onehalf[i-1] - K_plus_onehalf[i] - soil.sink[i]
end

i = soil.nsoi
a[i] = -K_plus_onehalf[i-1] / (2.0 * soil.dz_plus_onehalf[i-1])
c[i] = 0.0
b[i] = soil.cap[i] * soil.dz[i] / dt - a[i] - c[i]
d[i] = soil.cap[i] * soil.dz[i] / dt * soil.psi[i] - a[i] * (soil.psi[i-1] - soil.psi[i]) + K_plus_onehalf[i-1] - soil.K[i] - soil.sink[i]

# Solve for psi at n+1


soil.psi = tridiagonal_solver(a, b, c, d, soil.nsoi)
# --- Check water balance()

soil.Q0 = -K_onehalf / (2.0 * dz_onehalf) * ((soil.psi0 - psi_n[1]) + (soil.psi0 - soil.psi[1])) - K_onehalf
soil.QN = -soil.K[soil.nsoi]

soil.dtheta = 0.0
for i = 1:soil.nsoi
   soil.dtheta = soil.dtheta + (soil.theta[i] - theta_n[i]) * soil.dz[i]
end

soil.err = soil.dtheta - (soil.QN - soil.Q0) * dt

return soil

end # function

# Computational Grid
function compute_grid_settings(soil::soil_struct)
   # Set the Computational Grid for the Solver
   # Soil layer thickness (cm)
   for i = 1:soil.nsoi
      soil.dz[i] = 1.0
   end


   # Soil depth [cm] at i+1/2 interface between layers i & i+1 [negative distance from surface]

   soil.z_plus_onehalf[1] = -soil.dz[1]
   for i = 2:soil.nsoi
      soil.z_plus_onehalf[i] = soil.z_plus_onehalf[i-1] - soil.dz[i]
   end

   # Soil depth [cm] at center of layer i [negative distance from surface]

   soil.z[1]  = 0.5 * soil.z_plus_onehalf[1]
   soil.fz[1] = 1 - soil.bi^abs(soil.z[1])
   for i = 2:soil.nsoi
      soil.z[i] = 0.5 * (soil.z_plus_onehalf[i-1] + soil.z_plus_onehalf[i])

      # Assign the root fraction here
      soil.fz[i] = 1 - soil.bi^abs(soil.z[i]) - soil.fz[i-1] # if z is in cm
   end

   # Thickness between between z[i] & z[i+1]

   for i = 1:soil.nsoi-1
      soil.dz_plus_onehalf[i] = soil.z[i] - soil.z[i+1]
   end
   soil.dz_plus_onehalf[soil.nsoi] = 0.5 * soil.dz[soil.nsoi]



   return soil
end # Function


end # module
