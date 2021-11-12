###############################################################################
#
# Calculate diffusive thermo and SIF
#
###############################################################################
"""
    diffusive_S(τ_dd::Array{FT},
                ρ_dd::Array{FT},
                S⁻::Array{FT},
                S⁺::Array{FT},
                boundary_top::Array{FT},
                boundary_bottom::Array{FT},
                rsoil::Array{FT}
    ) where {FT<:AbstractFloat}

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
    # Get dimensions (1st is wavelength, 2nd is layers), for Stefab Boltzmann
    # just one effective wavelength
    nWL,nl = size(τ_dd);
    Xdd    = similar(τ_dd);
    Rdd    = zeros(FT,nWL,nl+1)
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
        net_diffuse[:,j]= (E⁻[:,j] + E⁺[:,j+1]).*(1 .-(τ_dd[:,j]+ρ_dd[:,j]))
    end

    return E⁻,E⁺,net_diffuse
end




"""
    diffusive_S!(
                sf_con::SFCache{FT},
                soil::SoilOpticals{FT},
                rt_dim::RTDimensions
    ) where {FT<:AbstractFloat}

Computes 2-stream diffusive radiation transport (used for thermal and SIF),
    given
- `sf_con` [`SFCache`](@ref) type cache
- `soil` [`SoilOpticals`](@ref) type struct
- `rt_dim` [`RTDimensions`](@ref) type struct
"""
function diffusive_S!(
            sf_con::SFCache{FT},
            soil::SoilOpticals{FT},
            rt_dim::RTDimensions
) where {FT<:AbstractFloat}
    # 1. unpack values from sf_con
    @unpack dnorm,  F⁻, F⁺, net_diffuse, Rdd, S⁻, S⁺, U, Xdd, Y, zeroB, ρ_dd,
            τ_dd = sf_con;
    @unpack ρ_SW_SIF = soil;
    @unpack nLayer, nLevel = rt_dim;

    # Get dimensions (1st is wavelength, 2nd is layers), for Stefan Boltzmann
    # just one effective wavelength
    dnorm .= 1 .- view(ρ_dd, :, 1);

    Rdd[:,end] .= ρ_SW_SIF;

    # Source at TOC (0 for SIF, downwelling LW from atmosphere for thermal)
    F⁻[:,1]  .= zeroB;
    # Source at bottom layer (0 for SIF, emitted from soil for thermal)
    U[:,end] .= zeroB;

    @inbounds for j=nLayer:-1:1   # bottom to top
        dnorm    .= 1 .- view(ρ_dd, :, j) .* view(Rdd, :, j+1);
        Xdd[:,j] .= view(τ_dd, :, j) ./ dnorm;
        Rdd[:,j] .= view(ρ_dd, :, j) .+ view(τ_dd, :, j) .* view(Rdd, :, j+1) .* view(Xdd, :, j);
        Y[:,j]   .= (view(ρ_dd, :, j) .* view(U, :, j+1) .+ view(S⁻, :, j)) ./ dnorm;
        U[:,j]   .= view(τ_dd, :, j) .* (view(Rdd, :, j+1) .* view(Y, :, j) .+ view(U, :, j+1)) .+ view(S⁺, :, j);
    end
    @inbounds for j=1:nLayer      # from top to bottom
        F⁻[:,j+1] .= view(Xdd, :, j) .* view(F⁻, :, j) .+ view(Y, :, j);
        F⁺[:,j  ] .= view(Rdd, :, j) .* view(F⁻, :, j) .+ view(U, :, j);
    end
    F⁺[:,end] .= view(Rdd, :, nLevel) .* view(F⁻, :, nLevel) .+ view(U, :, nLevel);

    # compute net diffuse radiation per layer:
    @inbounds for j=1:nLayer
        net_diffuse[:,j] .= (view(F⁻, :, j) .+ view(F⁺, :, j+1)) .* (1 .- (view(τ_dd, :, j) .+ view(ρ_dd, :, j)));
    end

    return nothing
end
