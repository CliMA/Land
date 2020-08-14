###############################################################################
#
# Simulate short wave radiation
#
###############################################################################
"""
    short_wave!(can::Canopy4RT{FT}, can_opt::CanopyOpticals, can_rad::CanopyRads{FT}, in_rad::IncomingRadiation{FT}, soil_opt::SoilOpticals{FT}, rt_con::RTContainer{FT}) where {FT<:AbstractFloat}

Simulate the short wave radiation through the canopy, given
- `can` [`Canopy4RT`](@ref) type struct of canopy information
- `can_opt` [`CanopyOpticals`](@ref) struct of optical layer properties
- `can_rad` A [`CanopyRads`](@ref) struct
- `in_rad` An [`IncomingRadiation`](@ref) struct
- `soil_opt` A [`SoilOpticals`](@ref) type struct for soil optical properties
- `rt_con` [`RTContainer`](@ref) type container
"""
function short_wave!(
            can::Canopy4RT{FT},
            can_opt::CanopyOpticals{FT},
            can_rad::CanopyRads{FT},
            in_rad::IncomingRadiation{FT},
            soil_opt::SoilOpticals{FT},
            rt_con::RTContainer{FT}
) where {FT<:AbstractFloat}
    # 1. unpack values from can and soil_opt
    @unpack LAI, nLayer, Ω = can;
    @unpack ks, sb, sf, sigb = can_opt;
    @unpack albedo_SW = soil_opt;
    sw_con = rt_con.sw_con;

    # 2. scattering and extinction coefficients to
    #    thin layer reflectances and transmittances
    # Eq. 17 in mSCOPE paper (changed here to compute real transmission)
    # TODO change to LAI distribution later
    iLAI = LAI * Ω / nLayer;
    τ_ss = exp(-ks * iLAI);
    sw_con.τ_dd .= 1 .- can_opt.a .* iLAI;
    sw_con.τ_sd .= sf   .* iLAI;
    sw_con.ρ_dd .= sigb .* iLAI;
    sw_con.ρ_sd .= sb   .* iLAI;
    @unpack ρ_dd, ρ_sd, τ_dd, τ_sd = sw_con;

    # 3. reflectance calculation
    # 3.1 Eq. 18 in mSCOPE paper
    Xss = τ_ss;

    # 3.2 Soil reflectance boundary condition (same for diffuse and direct)
    can_opt.R_sd[:,end] .= albedo_SW;
    can_opt.R_dd[:,end] .= albedo_SW;

    # 3.3 reflectance for each layer from bottom to top
    @inbounds for j in nLayer:-1:1
        sw_con.dnorm      .= 1 .-
                             view(ρ_dd, :, j) .*
                             view(can_opt.R_dd, :, j+1);
        can_opt.Xsd[:,j]  .= ( view(τ_sd, :, j) .+
                               Xss .* view(can_opt.R_sd, :, j+1) .*
                                      view(ρ_dd, :, j) ) ./
                             sw_con.dnorm;
        can_opt.Xdd[:,j]  .= view(τ_dd, :, j) ./ sw_con.dnorm;
        can_opt.R_sd[:,j] .= view(ρ_sd, :, j) .+
                             view(τ_dd, :, j) .*
                                ( view(can_opt.R_sd, :, j+1) .*
                                  Xss .+
                                  view(can_opt.R_dd, :, j+1) .*
                                  view(can_opt.Xsd , :, j  ) );
        can_opt.R_dd[:,j] .= view(ρ_dd        , :, j  ) .+
                             view(τ_dd        , :, j  ) .*
                             view(can_opt.R_dd, :, j+1) .*
                             view(can_opt.Xdd , :, j  );
    end

    # 4. flux profile calculation
    # Eq. 19 in mSCOPE paper
    # 4.1 Boundary condition at top: Incoming solar radiation
    can_opt.Es_[:,1]    .= in_rad.E_direct;
    can_rad.E_down[:,1] .= in_rad.E_diffuse;

    # 4.2 from top to bottom
    @inbounds for j=1:nLayer
        can_rad.netSW_sunlit[:,j] .= view(can_opt.Es_   , :, j) .*
                                     ( 1 .- (τ_ss .+ view(τ_sd, :, j) .+
                                     view(ρ_sd          , :, j)) );
        can_opt.Es_[:,j+1]        .= Xss .*
                                     view(can_opt.Es_   , :, j);
        can_rad.E_down[:,j+1]     .= view(can_opt.Xsd   , :, j) .*
                                     view(can_opt.Es_   , :, j) .+
                                     view(can_opt.Xdd   , :, j) .*
                                     view(can_rad.E_down, :, j);
        can_rad.E_up[:,j]         .= view(can_opt.R_sd  , :, j) .*
                                     view(can_opt.Es_   , :, j) .+
                                     view(can_opt.R_dd  , :, j) .*
                                     view(can_rad.E_down, :, j);
    end

    # 4.3 Boundary condition at the bottom, soil reflectance (Lambertian here)
    last_ind_co = lastindex(can_opt.R_sd, 2);
    can_rad.E_up[:,end] .= view(can_opt.R_sd  , :, last_ind_co) .*
                           view(can_opt.Es_   , :, last_ind_co) .+
                           view(can_opt.R_dd  , :, last_ind_co) .*
                           view(can_rad.E_down, :, last_ind_co);

    # 4.4 Hemispheric total outgoing
    can_rad.Eout .= view(can_rad.E_up, :, 1);

    # 4.5 compute net diffuse radiation per layer:
    @inbounds for j in 1:nLayer
        can_rad.netSW_shade[:,j] .= ( view(can_rad.E_down, :, j) .+
                                      view(can_rad.E_up, :, j+1) ) .*
                                    ( 1 .- ( view(τ_dd, :, j) .+
                                             view(ρ_dd, :, j) ) );
        # Add diffuse radiation to direct radiation as well:
        #can_rad.netSW_sunlit[:,j] += can_rad.netSW_shade[:,j]
    end

    # 4.6 outgoing in viewing direction
    # From Canopy
    sw_con.piLoc2 .= can_opt.vb .* view(can_opt.Po       , 1:nLayer)' .*
                                   view(can_rad.E_down, :, 1:nLayer)  .+
                     can_opt.vf .* view(can_opt.Po       , 1:nLayer)' .*
                                   view(can_rad.E_up  , :, 1:nLayer)  .+
                     can_opt.w  .* view(can_opt.Pso      , 1:nLayer)' .*
                                   in_rad.E_direct;
    #sw_con.piLoc  .= iLAI .* view(sum(sw_con.piLoc2, dims=2), :, 1);
    @inbounds for j in eachindex(sw_con.piLoc)
        sw_con.piLoc[j] = iLAI * sum( view(sw_con.piLoc2, j, :) );
    end

    # 4.7 From Soil
    sw_con.piLos .= view(can_rad.E_up, :, last_ind_co) .*
                    can_opt.Po[end];
    sw_con.piLo  .= sw_con.piLoc .+ sw_con.piLos;
    can_rad.Lo   .= sw_con.piLo ./ pi;

    # 4.8 Save albedos (hemispheric direct and diffuse and directional (obs))
    # rso and rdo are not computed separately
    can_rad.alb_obs     .= sw_con.piLo ./ ( in_rad.E_direct .+
                                            in_rad.E_diffuse );
    can_rad.alb_direct  .= view(can_opt.R_sd, :, 1);
    can_rad.alb_diffuse .= view(can_opt.R_dd, :, 1);

    return nothing
end
