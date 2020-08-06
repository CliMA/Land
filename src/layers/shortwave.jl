###############################################################################
#
# Simulate short wave radiation
#
###############################################################################
"""
    short_wave!(can::Canopy4RT{FT}, can_opt::CanopyOpticals, can_rad::CanopyRads{FT}, in_rad::IncomingRadiation{FT}, soil_opt::SoilOpticals{FT}) where {FT<:AbstractFloat}

Simulate the short wave radiation through the canopy, given
- `can` A [`Canopy4RT`](@ref) type struct for providing LAI and nLayer and clumping
- `can_opt` A [`CanopyOpticals`](@ref) struct for providing optical layer properties
- `can_rad` A [`CanopyRads`](@ref) struct
- `in_rad` An [`IncomingRadiation`](@ref) struct
- `soil_opt` A [`SoilOpticals`](@ref) type struct for soil optical properties
"""
function short_wave!(
            can::Canopy4RT{FT},
            can_opt::CanopyOpticals{FT},
            can_rad::CanopyRads{FT},
            in_rad::IncomingRadiation{FT},
            soil_opt::SoilOpticals{FT}
) where {FT<:AbstractFloat}
    # TODO Organize the code for better understanding
    # Unpack variables from can structure
    @unpack LAI, nLayer, Ω = can

    # TODO change to LAI distribution later
    iLAI = LAI * Ω / nLayer

    # Define soil as polynomial (depends on state vector size), still TBD in the structure mode now.
    # TODO emissivity of soil (goes into structure later)
    rsoil           = soil_opt.albedo_SW
    soil_emissivity = 1 .- rsoil

    # 2.2. convert scattering and extinction coefficients into thin layer reflectances and transmittances
    # Eq. 17 in mSCOPE paper (changed here to compute real transmission)
    τ_ss = exp(-can_opt.ks * iLAI)
    # Here, transmission looks ok, as it has to be sigf for iLAI=1!
    τ_dd = 1 .- can_opt.a * iLAI
    τ_sd = can_opt.sf   * iLAI
    ρ_sd = can_opt.sb   * iLAI
    ρ_dd = can_opt.sigb * iLAI

    # 2.3. reflectance calculation
    # Eq. 18 in mSCOPE paper
    # Just a scalar here but better to write it out:
    Xss  = τ_ss

    # Soil reflectance boundary condition (same for diffuse and direct)
    can_opt.R_sd[:,end] = rsoil
    can_opt.R_dd[:,end] = rsoil

    # Eq. 18 in mSCOPE paper
    # from bottom to top:
    @inbounds for j=nLayer:-1:1
        dnorm        = 1 .- ρ_dd[:,j] .* can_opt.R_dd[:,j+1]
        can_opt.Xsd[:,j]  = ( τ_sd[:,j] + Xss * can_opt.R_sd[:,j+1] .* ρ_dd[:,j] ) ./ dnorm
        can_opt.Xdd[:,j]  = τ_dd[:,j] ./ dnorm
        can_opt.R_sd[:,j] = ρ_sd[:,j] + τ_dd[:,j] .* ( can_opt.R_sd[:,j+1] * Xss + can_opt.R_dd[:,j+1] .* can_opt.Xsd[:,j] )
        can_opt.R_dd[:,j] = ρ_dd[:,j] + τ_dd[:,j] .* can_opt.R_dd[:,j+1] .* can_opt.Xdd[:,j]
    end

    # 3.2 flux profile calculation
    # Eq. 19 in mSCOPE paper
    # Boundary condition at top: Incoming solar radiation
    #in_rad.E_diffuse = in_rad.E_diffuse.+0.2# Add some here for checking, is almost 0 otherwise!
    can_opt.Es_[:,1]    = in_rad.E_direct
    can_rad.E_down[:,1] = in_rad.E_diffuse

    @inbounds for j=1:nLayer # from top to bottom
        can_rad.netSW_sunlit[:,j] = can_opt.Es_[:,j] .* ( 1 .- (τ_ss .+ τ_sd[:,j] + ρ_sd[:,j]) )
        can_opt.Es_[:,j+1]        = Xss .* can_opt.Es_[:,j]
        can_rad.E_down[:,j+1]     = can_opt.Xsd[:,j]  .* can_opt.Es_[:,j] + can_opt.Xdd[:,j]  .* can_rad.E_down[:,j]
        can_rad.E_up[:,j]         = can_opt.R_sd[:,j] .* can_opt.Es_[:,j] + can_opt.R_dd[:,j] .* can_rad.E_down[:,j]
    end
    # Boundary condition at the bottom, soil reflectance (Lambertian here)
    can_rad.E_up[:,end] = can_opt.R_sd[:,end] .* can_opt.Es_[:,end] + can_opt.R_dd[:,end] .* can_rad.E_down[:,end]

    # Hemispheric total outgoing (needs to go into stucture later):
    can_rad.Eout[:] = can_rad.E_up[:,1]

    # compute net diffuse radiation per layer:
    @inbounds for j=1:nLayer
        E_ = can_rad.E_down[:,j] + can_rad.E_up[:,j+1]
        can_rad.netSW_shade[:,j] = E_ .* ( 1 .- (τ_dd[:,j] + ρ_dd[:,j]) )
        # Add diffuse radiation to direct radiation as well:
        #can_rad.netSW_sunlit[:,j] += can_rad.netSW_shade[:,j]
    end

    # outgoing in viewing direction
    # From Canopy
    piLoc_  = iLAI * sum( can_opt.vb .* can_opt.Po[1:nLayer]' .* can_rad.E_down[:,1:nLayer] + can_opt.vf .* can_opt.Po[1:nLayer]' .* can_rad.E_up[:,1:nLayer] + can_opt.w .* can_opt.Pso[1:nLayer]' .* in_rad.E_direct, dims=2 )[:,1]
    # From Soil
    piLos_  = can_rad.E_up[:,end] * can_opt.Po[end]
    piLo_   = piLoc_ + piLos_
    can_rad.Lo   = piLo_  / pi
    # Save albedos (hemispheric direct and diffuse as well as directional (obs))
    can_rad.alb_obs     = piLo_ ./ ( in_rad.E_direct + in_rad.E_diffuse ) # rso and rdo are not computed separately
    can_rad.alb_direct  = can_opt.R_sd[:,1]
    can_rad.alb_diffuse = can_opt.R_dd[:,1]

    return nothing
end
