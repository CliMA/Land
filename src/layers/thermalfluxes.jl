###############################################################################
#
# Simulate thermal fluxes
#
###############################################################################
"""
    thermal_fluxes!(leaf_array::Array{LeafBios{FT},1}, can_opt::CanopyOpticals{FT}, can_rad::CanopyRads{FT}, can::Canopy4RT{FT}, soil_opt::SoilOpticals{FT}, incLW::Array{FT}, wl_set::WaveLengths{FT}) where {FT<:AbstractFloat}

Computes 2-stream diffusive radiation transport for thermal radiation (calls `compute_diffusive_S` internally).
Layer reflectance and transmission is computed from LW optical properties, layer sources from temperature and Planck law, boundary conditions from the atmosphere and soil emissivity and temperature.
Currently only uses Stefan Boltzmann law to compute spectrally integrated LW but can be easily adjusted to be spectrally resolved.
- `leaf_array` An array of [`LeafBios`](@ref) type struct (i.e. leaf optical properties can change with canopy height)
- `can_opt` A [`CanopyOpticals`](@ref) struct for providing optical layer properties
- `can_rad` A [`CanopyRads`](@ref) struct
- `can` A [`Canopy4RT`](@ref) type struct for providing LAI and nLayer and clumping
- `sO` A [`SoilOpticals`](@ref) type struct for soil optical properties
- `incLW` A 1D array with incoming long-wave radiation
- `wl_set` An [`WaveLengths`](@ref) type struct
"""
function thermal_fluxes!(
            leaf_array::Array{LeafBios{FT},1},
            can_opt::CanopyOpticals{FT},
            can_rad::CanopyRads{FT},
            can::Canopy4RT{FT},
            soil_opt::SoilOpticals{FT},
            incLW::Array{FT},
            wl_set::WaveLengths{FT}
) where {FT<:AbstractFloat}
    @unpack Ps, Po, Pso,ddf,ddb = can_opt
    @unpack T_sun, T_shade = can_rad
    @unpack Ω,nLayer,LAI,lidf = can
    @unpack albedo_LW, soil_skinT = soil_opt
    @unpack nwl = wl_set

    # Number of layers
    nl      = nLayer
    iLAI    = Ω*LAI/nLayer;

    # Sunlit fraction at mid layer:
    fSun    = (Ps[1:nLayer] + Ps[2:nLayer+1])/2;
    ϵ       = similar(T_shade)
    τ_dd    = zeros(FT,1,length(T_shade))
    ρ_dd    = similar(τ_dd)
    S⁺      = similar(ρ_dd)
    S⁻      = similar(ρ_dd)
    S_shade = similar(ρ_dd)
    S_sun   = similar(ρ_dd)

    # If we just have the same leaf everywhere, compute emissivities:
    if length(leaf_array)==1
        le   = leaf_array[1]
        # Compute layer properties:
        sigf = ddf*le.ρ_LW + ddb*le.τ_LW
        sigb = ddb*le.ρ_LW + ddf*le.τ_LW
        τ_dd = (1 - (1-sigf)*iLAI)*ones(nwl,nl)
        ρ_dd = (sigb*iLAI)*ones(nwl,nl)
        ϵ   .= (1 - τ_dd-ρ_dd);
    elseif length(leaf_array)==nLayer
        for i=1:nLayer
            le      = leaf_array[i]
            sigf    = ddf*le.ρ_LW + ddb*le.τ_LW
            sigb    = ddb*le.ρ_LW + ddf*le.τ_LW
            τ_dd[i] = (1 - (1-sigf)*iLAI)
            ρ_dd[i] = (sigb*iLAI)
            ϵ[i]    = (1 - τ_dd[i]-ρ_dd[i]);
        end
    else
        println("Complain, Array of leaves is neither 1 nor nLayer ")
    end

    # Only one wavelength --> do Stefan Boltzmann:
    #if length(wl)==1
    # Let's just do SB for now:
    if 1==1
        # Shaded leaves first, simple 1D array:
        S_shade= K_BOLTZMANN(FT) .* ϵ .* (T_shade.^4)
        # Sunlit leaves:
        if ndims(T_sun)>1
            @inbounds for i=1:length(T_shade)
                emi      = K_BOLTZMANN(FT) * ϵ[i] * T_sun[:,:,i].^4
                # weighted average over angular distribution
                S_sun[i] = mean(emi'*lidf);
            end
        else
            # Sunlit, simple 1D array:
            S_sun = K_BOLTZMANN(FT) .* ϵ .* T_sun.^4
        end
    else
        # Do Planck curve, tbd
    end
    S⁺[:] = iLAI*(fSun.*S_sun+(1 .-fSun).*S_shade)
    S⁻[:] = S⁺[:]
    soilEmission = K_BOLTZMANN(FT) * (1 .- albedo_LW) * soil_skinT^4
    # Run RT:
    F⁻,F⁺,net_diffuse = diffusive_S(τ_dd, ρ_dd,S⁻, S⁺,incLW, soilEmission, albedo_LW)
    for j = 1:nLayer
        can_rad.intNetLW_sunlit[j]  = net_diffuse[j]- 2S_sun[j];          #  sunlit leaf
        can_rad.intNetLW_shade[j]   = net_diffuse[j]- 2S_shade[j];     # shaded leaf
    end
    # Net soil LW as difference between up and downwelling at lowest level
    can_rad.RnSoilLW = F⁻[1,end]-F⁺[1,end]
    #@show F⁻[:,end]
    #@show F⁺[:,end]

    return F⁻,F⁺,net_diffuse
end
