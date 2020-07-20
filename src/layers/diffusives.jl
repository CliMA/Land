###############################################################################
#
# Calculate diffusive thermo and SIF
#
###############################################################################
"""
    diffusive_S(τ_dd::Array{FT}, ρ_dd::Array{FT}, S⁻::Array{FT}, S⁺::Array{FT}, boundary_top::Array{FT}, boundary_bottom::Array{FT}, rsoil::Array{FT}) where {FT<:AbstractFloat}

Computes 2-stream diffusive radiation transport (used for thermal and SIF) given:
- `τ_dd` A 2D Array with layer reflectances
- `ρ_dd` A 2D Array with layer transmissions
- `S⁻` A 2D Array with layer source terms in the downwelling direction
- `S⁺` A 2D Array with layer source terms in the upwelling direction
- `boundary_top` A 1D array with downwelling radiation at the top (top of canopy)
- `boundary_bottom` A 1D array with upwnwelling radiation at the bottom (soil)
- `rsoil` A 1D array with soil reflectance
"""
function diffusive_S(
            τ_dd::Array{FT},
            ρ_dd::Array{FT},
            S⁻::Array{FT},
            S⁺::Array{FT},
            boundary_top::Array{FT},
            boundary_bottom::Array{FT},
            rsoil::Array{FT}
) where {FT<:AbstractFloat}
    #Get dimensions (1st is wavelength, 2nd is layers), for Stefab Boltzmann, just one effective wavelength
    nwl,nl = size(τ_dd);
    Xdd    = similar(τ_dd);
    Rdd    = zeros(FT,nwl,nl+1)
    # Create Y,U matrices (see mSCOPE paper, eqs 30-39)
    Y  = similar(τ_dd);
    U  = similar(Rdd);
    E⁻ = similar(Rdd);
    E⁺ = similar(Rdd);
    net_diffuse = similar(τ_dd);

    Rdd[:,end] = rsoil

    # Source at TOC (0 for SIF, downwelling LW from atmosphere for thermal)
    E⁻[:,1]  = boundary_top;
    # Source at bottom layer (0 for SIF, emitted from soil for thermal)
    U[:,end].= boundary_bottom;

    @inbounds for j=nl:-1:1   # bottom to top
        dnorm    = 1 .-ρ_dd[:,j].*Rdd[:,j+1];
        Xdd[:,j] = τ_dd[:,j]./dnorm;
        Rdd[:,j] = ρ_dd[:,j]+τ_dd[:,j].*Rdd[:,j+1].*Xdd[:,j];
        Y[:,j]   = (ρ_dd[:,j].*U[:,j+1]+S⁻[:,j])./dnorm;
        U[:,j]   = τ_dd[:,j].*(Rdd[:,j+1].*Y[:,j]+U[:,j+1])+S⁺[:,j];
    end
    @inbounds for j=1:nl      # from top to bottom
        E⁻[:,j+1]= Xdd[:,j].*E⁻[:,j]+Y[:,j];
        E⁺[:,j]  = Rdd[:,j].*E⁻[:,j]+U[:,j];
    end
    E⁺[:,end] = Rdd[:,end].*E⁻[:,end]+U[:,end];
    # compute net diffuse radiation per layer:
    @inbounds for j=1:nl
        net_diffuse[:,j]= (E⁻[:,j] +  E⁺[:,j+1]).*(1 .-(τ_dd[:,j]+ρ_dd[:,j]))
    end

    return E⁻,E⁺,net_diffuse
end
