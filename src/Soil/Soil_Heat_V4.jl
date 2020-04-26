module Soil_Heat_V4
using Parameters

export physcon,
       soil,
       soil_thermal_properties,
       phase_change,
       soil_temperature,
       tridiagonal_solver,
       compute_grid_settings

# --- Physical constants in physcon structure
@with_kw mutable struct physcon{}
	tfrz = 273.15;                         # Freezing point of water [K]
	cwat = 4188.0;                         # Specific heat of water [J/kg/K]
	cice = 2117.27;                        # Specific heat of ice [J/kg/K]
	rhowat = 1000.0;                       # Density of water [kg/m3]
	rhoice = 917.0;                        # Density of ice [kg/m3]
	cvwat = cwat * rhowat; 					   # Heat capacity of water [J/m3/K]
	cvice = cice * rhoice; 					   # Heat capacity of ice [J/m3/K]
	tkwat = 0.57;                          # Thermal conductivity of water [W/m/K]
	tkice = 2.29;                          # Thermal conductivity of ice [W/m/K]
	hfus = 0.3337e6;                       # Heat of fusion for water at 0 C [J/kg]
end

# Soil texture classes [Cosby et al. 1984. Water Resources Research 20:682-690]

#  1: sand
#  2: loamy sand
#  3: sandy loam
#  4: silty loam
#  5: loam
#  6: sandy clay loam
#  7  silty clay loam
#  8: clay loam
#  9: sandy clay
# 10: silty clay
# 11: clay

@with_kw mutable struct soil{}
	silt = [ 5.0, 12.0, 32.0, 70.0, 39.0, 15.0, 56.0, 34.0,  6.0, 47.0, 20.0]; # Percent silt
	sand = [92.0, 82.0, 58.0, 17.0, 43.0, 58.0, 10.0, 32.0, 52.0,  6.0, 22.0]; # Percent sand
	clay = [ 3.0,  6.0, 10.0, 13.0, 18.0, 27.0, 34.0, 34.0, 42.0, 47.0, 58.0]; # Percent clay

	# Volumetric soil water content at saturation [porosity]
	# (Clapp & Hornberger. 1978. Water Resources Research 14:601-604)

	watsat = [0.395, 0.410, 0.435, 0.485, 0.451, 0.420, 0.477, 0.476, 0.426, 0.492, 0.482]

	nsoi            # Number of soil layers
	dz              = ones(nsoi)    # Soil layer thickness (cm)
	z_plus_onehalf  = zeros(nsoi) 	# Soil depth [cm] at i+1/2 interface between layers i & i+1 [negative distance from surface]
	z               = zeros(nsoi) 	# Soil depth [cm] at center of layer i [negative distance from surface]
	dz_plus_onehalf = zeros(nsoi) 	# Thickness between between z[i] & z[i+1]
	soil_texture    = 1 		  	# Initial Soil Texture Class
	method          = "excess-heat"	# Initial method for Phase Change
	tsoi            = zeros(nsoi) 	# Soil Temperature Grid
	gsoi 			= 0   			# Ground Heat Flux
	hfsoi 			= 0				# Initialize total soil heat of fusion to zero
	h2osoi_ice		= zeros(nsoi) 	# Actual water content as ice fraction
	h2osoi_liq		= zeros(nsoi) 	# Actual water content as liquid fraction
	tk 				= zeros(nsoi) 	# Thermal conductivity
	cv 				= zeros(nsoi) 	# Heat Capacity

	tsoi0 			= zeros(nsoi) 	# Saving current soil temperature for energy conservation check
	tk_plus_onehalf = zeros(nsoi) 	# Thermal conductivity at interface [W/m/K]

end

# Computational Grid
function compute_grid_settings(soil::soil)
   # Set the Computational Grid for the Solver
   # Soil layer thickness (m)
   for i = 1:soil.nsoi
      soil.dz[i] = 0.025
   end


   # Soil depth [cm] at i+1/2 interface between layers i & i+1 [negative distance from surface]

   soil.z_plus_onehalf[1] = -soil.dz[1]
   for i = 2:soil.nsoi
      soil.z_plus_onehalf[i] = soil.z_plus_onehalf[i-1] - soil.dz[i]
   end

   # Soil depth [cm] at center of layer i [negative distance from surface]

   soil.z[1]  = 0.5 * soil.z_plus_onehalf[1]
   for i = 2:soil.nsoi
      soil.z[i] = 0.5 * (soil.z_plus_onehalf[i-1] + soil.z_plus_onehalf[i])
   end

   # Thickness between between z[i] & z[i+1]

   for i = 1:soil.nsoi-1
      soil.dz_plus_onehalf[i] = soil.z[i] - soil.z[i+1]
   end
   soil.dz_plus_onehalf[soil.nsoi] = 0.5 * soil.dz[soil.nsoi]

   return soil
end # Function

# Function for soil thermal conductivity and specific heat
function soil_thermal_properties(physcon::physcon, soil::soil)

# Calculate soil thermal conductivity & heat capacity

# ------------------------------------------------------
# Input
#   physcon.hfus             ! Heat of fusion for water at 0 C [J/kg]
#   physcon.tfrz             ! Freezing point of water [K]
#   physcon.tkwat            ! Thermal conductivity of water [W/m/K]
#   physcon.tkice            ! Thermal conductivity of ice [W/m/K]
#   physcon.cvwat            ! Heat capacity of water [J/m3/K]
#   physcon.cvice            ! Heat capacity of ice [J/m3/K]
#   physcon.rhowat           ! Density of water [kg/m3]
#   physcon.rhoice           ! Density of ice [kg/m3]
#   soil.method           ! Use excess heat | apparent heat capacity for phase change
#   soil.soil_texture     ! Soil texture class()
#   soil.sand             ! Percent sand
#   soil.watsat           ! Volumetric soil water content at saturation [porosity]
#   soil.nsoi             ! Number of soil layers
#   soil.dz               ! Soil layer thickness [m]
#   soil.tsoi             ! Soil temperature [K]
#   soil.h2osoi_liq       ! Unfrozen water, liquid [kg H2O/m2]
#   soil.h2osoi_ice       ! Frozen water, ice [kg H2O/m2]
#
# Input/output
#   soil.tk               ! Thermal conducitivty [W/m/K]
#   soil.cv               ! Volumetric heat capacity [J/m3/K]
# ------------------------------------------------------

for i = 1:soil.nsoi

   # --- Soil texture to process

   k = soil.soil_texture

   # --- Volumetric soil water & ice

   watliq = soil.h2osoi_liq[i] / (physcon.rhowat * soil.dz[i])
   watice = soil.h2osoi_ice[i] / (physcon.rhoice * soil.dz[i])

   # Fraction of total volume that is liquid water

   fliq = watliq / (watliq + watice)

   # Soil water relative to saturation

   s = min((watliq + watice)/soil.watsat[k], 1)

   # --- Dry thermal conductivity [W/m/K] from bulk density [kg/m3]

   bd = 2700 * (1 - soil.watsat[k])
   tkdry = (0.135 * bd + 64.7) / (2700 - 0.947 * bd)

   # --- Soil solids thermal conducitivty [W/m/K]

   # Thermal conductivity of quartz [W/m/K]

   tk_quartz = 7.7

   # Quartz fraction

   quartz = soil.sand[k] / 100

   # Thermal conductivity of other minerals [W/m/K]

   if (quartz > 0.2)
      tko = 2
   else()
      tko = 3
   end

   # Thermal conductivity of soil solids [W/m/K]

   tksol = tk_quartz^quartz * tko^(1-quartz)

   # --- Saturated thermal conductivity [W/m/K] and unfrozen & frozen values

   tksat = tksol^(1-soil.watsat[k]) * physcon.tkwat^(fliq*soil.watsat[k]) * physcon.tkice^(soil.watsat[k]-fliq*soil.watsat[k])
   tksat_u = tksol^(1-soil.watsat[k]) * physcon.tkwat^soil.watsat[k]
   tksat_f = tksol^(1-soil.watsat[k]) * physcon.tkice^soil.watsat[k]

   # --- Kersten number and unfrozen & frozen values

   if (soil.sand[k] < 50)
      ke_u = log10(max(s,0.1)) + 1
   else()
      ke_u = 0.7 * log10(max(s,0.05)) + 1
   end
   ke_f = s

   if (soil.tsoi[i] >= physcon.tfrz)
      ke = ke_u
   else()
      ke = ke_f
   end
   # --- Thermal conductivity [W/m/K] and unfrozen & frozen values

   soil.tk[i] = (tksat - tkdry) * ke + tkdry
   tku = (tksat_u - tkdry) * ke_u + tkdry
   tkf = (tksat_f - tkdry) * ke_f + tkdry

   # --- Heat capacity of soil solids [J/m3/K]

   cvsol = 1.926e06

   # --- Heat capacity [J/m3/K] and unfrozen & frozen values

   soil.cv[i] = (1 - soil.watsat[k]) * cvsol + physcon.cvwat * watliq + physcon.cvice * watice
   cvu = (1 - soil.watsat[k]) * cvsol + physcon.cvwat * (watliq + watice)
   cvf = (1 - soil.watsat[k]) * cvsol + physcon.cvice * (watliq + watice)

   # --- Adjust heat capacity & thermal conductivity if using apparent heat capacity

   if soil.method == "apparent-heat-capacity"

      # Temperature range for freezing & thawing [K]

      tinc = 0.5

      # Heat of fusion [J/m3] - This is equivalent to ql = hfus * (h2osoi_liq + h2osoi_ice) / dz

      ql = physcon.hfus * (physcon.rhowat * watliq + physcon.rhoice * watice)

      # Heat capacity & thermal conductivity

      if (soil.tsoi[i] > physcon.tfrz+tinc)
         soil.cv[i] = cvu
         soil.tk[i] = tku
      end

      if (soil.tsoi[i] >= physcon.tfrz-tinc && soil.tsoi[i] <= physcon.tfrz+tinc)
         soil.cv[i] = (cvf + cvu) / 2 + ql / (2 * tinc)
         soil.tk[i] = tkf + (tku - tkf) * (soil.tsoi[i] - physcon.tfrz + tinc) / (2 * tinc)
      end

      if (soil.tsoi[i] < physcon.tfrz-tinc)
         soil.cv[i] = cvf
         soil.tk[i] = tkf
      end
   end

end

return soil
end # function

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

function phase_change(physcon::physcon, soil::soil, dt)

# Adjust temperatures for phase change. Freeze | melt ice using
# energy excess | deficit needed to change temperature to the
# freezing point.

# ------------------------------------------------------
# Input
#   dt                      ! Time step [s]
#   physcon.hfus            ! Heat of fusion for water at 0 C [J/kg]
#   physcon.tfrz            ! Freezing point of water [K]
#   soil.nsoi            ! Number of soil layers
#   soil.dz              ! Soil layer thickness [m]
#   soil.cv              ! Volumetric heat capacity [J/m3/K]
#
# Input/output
#   soil.tsoi            ! Soil temperature [K]
#   soil.h2osoi_liq      ! Unfrozen water, liquid [kg H2O/m2]
#   soil.h2osoi_ice      ! Frozen water, ice [kg H2O/m2]
#
# Output
#   soil.hfsoi           ! Soil phase change energy flux [W/m2]
# ------------------------------------------------------

# --- Initialize total soil heat of fusion to zero

soil.hfsoi = 0

# --- Now loop over all soil layers to calculate phase change

for i = 1:soil.nsoi

   # --- Save variables prior to phase change

   wliq0 = soil.h2osoi_liq[i];     # Amount of liquid water before phase change
   wice0 = soil.h2osoi_ice[i];     # Amount of ice before phase change
   wmass0 = wliq0 + wice0;            # Amount of total water before phase change
   tsoi0 = soil.tsoi[i];           # Soil temperature before phase change

   # --- Identify melting | freezing layers & set temperature to freezing

   # Default condition is no phase change [imelt = 0]

   imelt = 0

   # Melting: if ice exists above melt point; melt some to liquid.
   # Identify melting by imelt = 1

   if (soil.h2osoi_ice[i] > 0 && soil.tsoi[i] > physcon.tfrz)
      imelt = 1
      soil.tsoi[i] = physcon.tfrz
   end

   # Freezing: if liquid exists below melt point; freeze some to ice.
   # Identify freezing by imelt = 2

   if (soil.h2osoi_liq[i] > 0 && soil.tsoi[i] < physcon.tfrz)
      imelt = 2
      soil.tsoi[i] = physcon.tfrz
   end

   # --- Calculate energy for freezing | melting

   # The energy for freezing | melting [W/m2] is assessed from the energy
   # excess | deficit needed to change temperature to the freezing point.
   # This is a potential energy flux; because cannot melt more ice than is()
   # present | freeze more liquid water than is present.
   #
   # heat_flux_pot .> 0: freezing; heat_flux_pot .< 0: melting

   if (imelt > 0)
      heat_flux_pot = (soil.tsoi[i] - tsoi0) * soil.cv[i] * soil.dz[i] / dt
   else()
      heat_flux_pot = 0
   end

   # Maximum energy for melting | freezing [W/m2]

   if (imelt == 1)
      heat_flux_max = -soil.h2osoi_ice[i] * physcon.hfus / dt
   end

   if (imelt == 2)
      heat_flux_max = soil.h2osoi_liq[i] * physcon.hfus / dt
   end

   # --- Now freeze | melt ice

   if (imelt > 0)

      # Change in ice [kg H2O/m2/s]: freeze [+] | melt [-]

      ice_flux = heat_flux_pot / physcon.hfus

      # Update ice [kg H2O/m2]

      soil.h2osoi_ice[i] = wice0 + ice_flux * dt

      # Cannot melt more ice than is present

      soil.h2osoi_ice[i] = max(0, soil.h2osoi_ice[i])

      # Ice cannot exceed total water that is present

      soil.h2osoi_ice[i] = min(wmass0, soil.h2osoi_ice[i])

      # Update liquid water [kg H2O/m2] for change in ice

      soil.h2osoi_liq[i] = max(0, (wmass0-soil.h2osoi_ice[i]))

      # Actual energy flux from phase change [W/m2]. This is equal to
      # heat_flux_pot except if tried to melt too much ice.

      heat_flux = physcon.hfus * (soil.h2osoi_ice[i] - wice0) / dt

      # Sum energy flux from phase change [W/m2]

      soil.hfsoi = soil.hfsoi + heat_flux

      # Residual energy not used in phase change is added to soil temperature

      residual = heat_flux_pot - heat_flux
      soil.tsoi[i] = soil.tsoi[i] - residual * dt / (soil.cv[i] * soil.dz[i])

      # Error check: make sure actual phase change does not exceed permissible phase change

      if (abs(heat_flux) > abs(heat_flux_max))
         error("Soil temperature energy conservation error: phase change")
      end

      # Freezing: make sure actual phase change does not exceed permissible phase change
      # & that the change in ice does not exceed permissible change

      if (imelt == 2)

         # Energy flux [W/m2]

         constraint = min(heat_flux_pot, heat_flux_max)
         err = heat_flux - constraint
         if (abs(err) > 1e-03)
            error("Soil temperature energy conservation error: freezing energy flux")
         end

         # Change in ice [kg H2O/m2]

         err = (soil.h2osoi_ice[i] - wice0) - constraint / physcon.hfus * dt
         if (abs(err) > 1e-03)
            error("Soil temperature energy conservation error: freezing ice flux")
         end
      end

      # Thawing: make sure actual phase change does not exceed permissible phase change
      # & that the change in ice does not exceed permissible change

      if (imelt == 1)

         # Energy flux [W/m2]

         constraint = max(heat_flux_pot, heat_flux_max)
         err = heat_flux - constraint
         if (abs(err) > 1e-03)
            error("Soil temperature energy conservation error: thawing energy flux")
         end

         # Change in ice [kg H2O/m2]

         err = (soil.h2osoi_ice[i] - wice0) - constraint / physcon.hfus * dt
         if (abs(err) > 1e-03)
            error("Soil temperature energy conservation error: thawing ice flux")
         end
      end

   end

end

return soil

end # function

function soil_temperature(physcon::physcon, soil::soil, tsurf, dt)

# Use an implicit formulation with the surface boundary condition specified
# as the surface temperature to solve for soil temperatures at time n+1.
#
# Calculate soil temperatures as:
#
#      dT   d     dT
#   cv -- = -- (k --)
#      dt   dz    dz
#
# where: T = temperature [K]
#        t = time [s]
#        z = depth [m]
#        cv = volumetric heat capacity [J/m3/K]
#        k = thermal conductivity [W/m/K]
#
# Set up a tridiagonal system of equations to solve for T at time n+1;
# where the temperature equation for layer i is()
#
#   d_i = a_i [T_i-1] n+1 + b_i [T_i] n+1 + c_i [T_i+1] n+1
#
# For soil layers undergoing phase change, set T_i = Tf [freezing] & use
# excess energy to freeze | melt ice:
#
#   Hf_i = (Tf - [T_i] n+1) * cv_i * dz_i / dt
#
# During the phase change; the unfrozen & frozen soil water
# (h2osoi_liq, h2osoi_ice) are adjusted.
#
# Or alternatively; use the apparent heat capacity method to
# account for phase change. In this approach; h2osoi_liq
# & h2osoi_ice are not calculated.
#
# ------------------------------------------------------
# Input
#   tsurf                   ! Surface temperature [K]
#   dt                      ! Time step [s]
#   soil.method          ! Use excess heat | apparent heat capacity for phase change
#   soil.nsoi            ! Number of soil layers
#   soil.z               ! Soil depth [m]
#   soil.z_plus_onehalf  ! Soil depth [m] at i+1/2 interface between layers i & i+1
#   soil.dz              ! Soil layer thickness [m]
#   soil.dz_plus_onehalf ! Thickness [m] between between i & i+1
#   soil.tk              ! Thermal conductivity [W/m/K]
#   soil.cv              ! Heat capacity [J/m3/K]
#
# Input/output
#   soil.tsoi            ! Soil temperature [K]
#   soil.h2osoi_liq      ! Unfrozen water, liquid [kg H2O/m2]
#   soil.h2osoi_ice      ! Frozen water, ice [kg H2O/m2]
#
# Output
#   soil.gsoi            ! Energy flux into soil [W/m2]
#   soil.hfsoi           ! Soil phase change energy flux [W/m2]
# ------------------------------------------------------

# --- Save current soil temperature for energy conservation check

for i = 1:soil.nsoi
   soil.tsoi0[i] = soil.tsoi[i]
end

# --- Thermal conductivity at interface [W/m/K]

for i = 1:soil.nsoi-1
   soil.tk_plus_onehalf[i] = soil.tk[i] * soil.tk[i+1] * (soil.z[i]-soil.z[i+1]) / (soil.tk[i]*(soil.z_plus_onehalf[i]-soil.z[i+1]) + soil.tk[i+1]*(soil.z[i]-soil.z_plus_onehalf[i]))
end

# --- Set up tridiagonal matrix

# Terms for tridiagonal matrix
a = zeros(soil.nsoi);
b = zeros(soil.nsoi);
c = zeros(soil.nsoi);
d = zeros(soil.nsoi);


# Top soil layer with tsurf as boundary condition

i = 1
m = soil.cv[i] * soil.dz[i] / dt
a[i] = 0
c[i] = -soil.tk_plus_onehalf[i] / soil.dz_plus_onehalf[i]
b[i] = m - c[i] + soil.tk[i] / (0 - soil.z[i])
d[i] = m * soil.tsoi[i] + soil.tk[i] / (0 - soil.z[i]) * tsurf

# Layers 2 to nsoi-1

for i = 2:soil.nsoi-1
   m = soil.cv[i] * soil.dz[i] / dt
   a[i] = -soil.tk_plus_onehalf[i-1] / soil.dz_plus_onehalf[i-1]
   c[i] = -soil.tk_plus_onehalf[i] / soil.dz_plus_onehalf[i]
   b[i] = m - a[i] - c[i]
   d[i] = m * soil.tsoi[i]
end

# Bottom soil layer with zero heat flux

i = soil.nsoi
m = soil.cv[i] * soil.dz[i] / dt
a[i] = -soil.tk_plus_onehalf[i-1] / soil.dz_plus_onehalf[i-1]
c[i] = 0
b[i] = m - a[i]
d[i] = m * soil.tsoi[i]

# --- Solve for soil temperature

soil.tsoi = tridiagonal_solver(a, b, c, d, soil.nsoi)

# --- Derive energy flux into soil [W/m2]

soil.gsoi = soil.tk[1] * (tsurf - soil.tsoi[1]) / (0 - soil.z[1])

# --- Phase change for soil layers undergoing freezing of thawing

if soil.method == "apparent-heat-capacity"

   # No explicit phase change energy flux. This is included in the heat capacity.

   soil.hfsoi = 0

   elseif  soil.method == "excess-heat"

   # Adjust temperatures for phase change. Freeze | melt ice using energy
   # excess | deficit needed to change temperature to the freezing point.
   # The variable hfsoi is returned as the energy flux from phase change [W/m2].

   soil = phase_change(physcon, soil, dt)

end

# --- Check for energy conservation

# Sum change in energy [W/m2]

edif = 0
for i = 1:soil.nsoi
   edif = edif + soil.cv[i] * soil.dz[i] * (soil.tsoi[i] - soil.tsoi0[i]) / dt
end

# Error check

err = edif - soil.gsoi - soil.hfsoi
if (abs(err) > 1e-03)
   error("Soil temperature energy conservation error")
end

return soil

end # function

end # module