###############################################################################
#
# Volume scattering functions adatpted from Volscatt version 2 by W. Verhoef
#
###############################################################################
"""
    volscatt!(cache::Array{FT,1},
              sza::FT,
              vza::FT,
              raa::FT,
              ttl::FT
    ) where {FT<:AbstractFloat}

Calculate interception parameters (`chi_s` and `chi_s`) and leaf reflectance
    multiplier (`frho`) and transmittance multiplier (`ftau`), given
- `cache` Array cache for results
- `sza` Solar zenith angle
- `vza` Viewing zenith angle
- `raa` Relative azimuth angle
- `ttl` Leaf inclination angle
"""
function volscatt!(
            cache::Array{FT,1},
            sza::FT,
            vza::FT,
            raa::FT,
            ttl::FT
) where {FT<:AbstractFloat}
    psi_rad = deg2rad(raa);
    cos_raa = cosd(raa);
    cos_ttl = cosd(ttl);
    sin_ttl = sind(ttl);
    cos_sza = cosd(sza);
    sin_sza = sind(sza);
    cos_vza = cosd(vza);
    sin_vza = sind(vza);
    Cs      = cos_ttl*cos_sza;
    Ss      = sin_ttl*sin_sza;
    Co      = cos_ttl*cos_vza;
    So      = sin_ttl*sin_vza;

    cosbts = FT(1.0);
    cosbto = FT(1.0);
    if (abs(Ss)>1e-6)
        cosbts = -Cs / Ss;
    end
    if (abs(So)>1e-6)
        cosbto = -Co / So;
    end

    if (abs(cosbts)<1)
        bts = acos(cosbts);
        ds  = Ss;
    else
        bts = FT(pi);
        ds  = Cs;
    end

    if abs(cosbto)<1
        bto = acos(cosbto);
        doo = So;
    elseif vza<90
        bto = FT(pi);
        doo = Co;
    else
        bto = 0;
        doo = -Co;
    end

    chi_s = 2 / FT(pi) * ( (bts - FT(pi)/2) * Cs + sin(bts) * Ss );
    chi_o = 2 / FT(pi) * ( (bto - FT(pi)/2) * Co + sin(bto) * So );

    # Computation of auxiliary azimut angles bt1, bt2, bt3 used
    # for the computation of the bidirectional scattering coefficient w

    btran1 = abs(bts - bto);
    btran2 = 2FT(pi) - bts - bto;
    #btran_=pi-abs(bts+bto-FT(pi))

    if psi_rad <= btran1
        bt1 = psi_rad;
        bt2 = btran1;
        bt3 = btran2;
    else
        bt1 = btran1;
        if psi_rad <= btran2
            bt2 = psi_rad;
            bt3 = btran2;
        else
            bt2 = btran2;
            bt3 = psi_rad;
        end
    end

    t1 = 2Cs * Co + Ss * So * cos_raa;
    t2 = 0;
    if bt2 > 0
        t2 = sin(bt2) * ( 2ds * doo + Ss * So * cos(bt1) * cos(bt3) );
    end

    denom = 2 * FT(pi)^2;
    frho  = ((pi-bt2) * t1 + t2) / denom;
    ftau  = (-bt2     * t1 + t2) / denom;

    # fill the values to cache
    cache[1] = chi_s;
    cache[2] = abs(chi_o);
    cache[3] = frho;
    cache[4] = ftau;

    return nothing
end
