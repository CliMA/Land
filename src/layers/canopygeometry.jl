###############################################################################
#
# Update clumping factor
#
###############################################################################
"""
    clumping_factor!(can::Canopy4RT{FT}, angles::SolarAngles{FT}) where {FT<:AbstractFloat}

Calculate the clumping factor, given
- `can` [`Canopy4RT`](@ref) type struct, providing canopy structure information
- `angles` [`SolarAngles`](@ref) type struct, defining sun-sensor geometry
"""
function clumping_factor!(
            can::Canopy4RT{FT},
            angles::SolarAngles{FT}
) where {FT<:AbstractFloat}
    @unpack clump_a, clump_b = can;
    @unpack tts = angles;

    if clump_b > 0
        can.Ω = clump_a + clump_b * (1 - cosd(tts));
    end

    return nothing
end








###############################################################################
#
# Update canopy geometry
#
###############################################################################
"""
    canopy_geometry!(can::Canopy4RT{FT}, angles::SolarAngles{FT}, can_opt::CanopyOpticals{FT}, rt_con::RTContainer{FT}) where {FT<:AbstractFloat}

Computes canopy optical properties (extinction coefficients for direct and
    diffuse light) based on the SAIL model. Most important input parameters are
    leaf inclination and azimuth distribution functions and sun-sensor
    geometry . Canopy clumping Ω is implemented as in Pinty et al (2015), given
- `can` [`Canopy4RT`](@ref) type struct, providing canopy structure information
- `angles` [`SolarAngles`](@ref) type struct, defining sun-sensor geometry
- `can_opt` [`CanopyOpticals`](@ref) type, where optical properties is stored
- `rt_con` [`RTContainer`](@ref) type container
"""
function canopy_geometry!(
            can::Canopy4RT{FT},
            angles::SolarAngles{FT},
            can_opt::CanopyOpticals{FT},
            rt_con::RTContainer{FT}
) where {FT<:AbstractFloat}
    # 1. update clumping factor from zenith angle
    clumping_factor!(can, angles);

    # 2. update solor angle dependent variables
    @unpack tts, tto, psi = angles;
    cos_tto = cosd(tto);
    tan_tto = tand(tto);
    cos_psi = cosd(psi);
    sin_tto = sind(tto);
    cos_tts = cosd(tts);
    sin_tts = sind(tts);
    tan_tts	= tand(tts);
    cts_cto = cos_tts * cos_tto;
    dso		= sqrt( tan_tts^2 + tan_tto^2 - 2*tan_tts*tan_tto*cos_psi );
    psi_vol = abs( psi - 360*round(psi/360) );

    # 3. unpack canopy parameters
    @unpack dx, hot, LAI, lazitab, lidf, litab, nLayer, xl, xl_e, Ω = can;

    # 4. update the RTContainer
    rt_con.cos_philo .= cosd.(lazitab .- psi);
    @unpack cos_philo, cos_ttli, cos_ttlo, sin_ttli = rt_con;

    # 5. calculate geometric factors associated with extinction and scattering
    can_opt.ks  = 0;
    can_opt.ko  = 0;
    can_opt.bf  = 0;
    can_opt.sob = 0;
    can_opt.sof = 0;
    @inbounds for i=1:length(litab)
        _lit = litab[i];
        _lid = lidf[i];
        _ctl = cos_ttli[i];

        # interception parameters and ref/trans multipliers
        volscatt!(rt_con.vol_scatt, tts, tto, psi_vol, _lit);
        chi_s, chi_o, frho, ftau = rt_con.vol_scatt;

        # Extinction coefficients
        ksli = abs(chi_s / cos_tts);
        koli = abs(chi_o / cos_tto);

        # Area scattering coefficient fractions
        if ftau < 0
            sobli = abs(ftau) * FT(pi) / cts_cto;
            sofli = abs(frho) * FT(pi) / cts_cto;
        else
            sobli = frho * FT(pi) / cts_cto;
            sofli = ftau * FT(pi) / cts_cto;
        end

        # add up the values in each layer
        can_opt.ks  += ksli   * _lid;
        can_opt.ko  += koli   * _lid;
        can_opt.bf  += _ctl^2 * _lid;
        can_opt.sob += sobli  * _lid;
        can_opt.sof += sofli  * _lid;
    end

    # 6. geometric factors to be used later with rho and tau
    @unpack bf, ko, ks = can_opt;
    can_opt.sdb = (ks + bf) / 2;
    can_opt.sdf = (ks - bf) / 2;
    can_opt.dob = (ko + bf) / 2;
    can_opt.dof = (ko - bf) / 2;
    can_opt.ddb = (1  + bf) / 2;
    can_opt.ddf = (1  - bf) / 2;

    # 7. eq 19 in vdT 2009 page 305 modified by Joris
    rt_con._Cs .= cos_ttli .* cos_tts; # [nli]
    rt_con._Ss .= sin_ttli .* sin_tts; # [nli]
    rt_con._Co .= cos_ttli .* cos_tto; # [nli]
    rt_con._So .= sin_ttli .* sin_tto; # [nli]
    @unpack _Co, _Cs, _So, _Ss, _1s = rt_con;
    # rt_con.cds .= _Cs * _1s .+ _Ss * cos_ttlo' ; # [nli, nlazi]
    # rt_con.cdo .= _Co * _1s .+ _So * cos_philo'; # [nli, nlazi]
    mul!(rt_con.cds, _Cs, _1s);
    mul!(rt_con._2d, _Ss, cos_ttlo' );
    rt_con.cds .+= rt_con._2d;
    mul!(rt_con.cdo, _Co, _1s);
    mul!(rt_con._2d, _So, cos_philo');
    rt_con.cdo .+= rt_con._2d;
    @unpack cdo, cds = rt_con;

    # 8. update fs and fo
    # This is basically equivalent to Kb in Bonan, eq. 14.21
    # TOD reduce allocations
    can_opt.fs      .= cds ./ cos_tts;
    can_opt.fo      .= cdo ./ cos_tto;
    can_opt.absfs   .= abs.( can_opt.fs );
    can_opt.cosΘ_l  .= cos_ttli .* _1s;
    can_opt.cos2Θ_l .= can_opt.cosΘ_l .^ 2;
    can_opt.fsfo    .= can_opt.fs .* can_opt.fo;
    can_opt.absfsfo .= abs.( can_opt.fsfo );

    # 9. probabilities Ps, Po, Pso
    _fac_s = exp(ks*Ω*LAI) * (1 - exp(-ks*Ω*LAI*dx)) / (ks*Ω*LAI*dx);
    _fac_o = exp(ko*Ω*LAI) * (1 - exp(-ko*Ω*LAI*dx)) / (ko*Ω*LAI*dx);
    can_opt.Ps  .= xl_e .* _fac_s;
    can_opt.Po  .= xl_e .* _fac_o;
    @inline f(x) = psofunction(ko, ks, Ω, LAI, hot, dso, x);

    # TODO minimize the allocations here
    # length(xl) * 7 allocations here!
    @inbounds for j=1:length(xl)
        can_opt.Pso[j] = quadgk(f, xl[j]-dx, xl[j], rtol=1e-2)[1] / dx;
    end

    # takes care of rounding error
    # can_opt.Pso[can_opt.Pso.>can_opt.Po] =
    # minimum([can_opt.Po[can_opt.Pso.>can_opt.Po]
    #          can_opt.Ps[can_opt.Pso.>can_opt.Po]],dims=2)
    # can_opt.Pso[can_opt.Pso.>can_opt.Ps] =
    # minimum([can_opt.Po[can_opt.Pso.>can_opt.Ps]
    #          can_opt.Ps[can_opt.Pso.>can_opt.Ps]],dims=2)

    return nothing
end
