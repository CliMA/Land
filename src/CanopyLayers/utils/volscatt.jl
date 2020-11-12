###############################################################################
#
# Volume scattering functions adatpted from Volscatt version 2 by W. Verhoef
#
###############################################################################
"""
    volscatt!(container::Array{FT,1},
              tts::FT,
              tto::FT,
              psi::FT,
              ttl::FT
    ) where {FT<:AbstractFloat}

Calculate interception parameters (`chi_s` and `chi_s`) and leaf reflectance
    multiplier (`frho`) and transmittance multiplier (`ftau`), given
- `container` Array container for results
- `tts` Solar zenith angle
- `tto` Viewing zenith angle
- `psi` Azimuth angle
- `ttl` Leaf inclination angle
"""
function volscatt!(
            container::Array{FT,1},
            tts::FT,
            tto::FT,
            psi::FT,
            ttl::FT
) where {FT<:AbstractFloat}
    psi_rad = deg2rad(psi);
    cos_psi = cosd(psi);
    cos_ttl = cosd(ttl);
    sin_ttl = sind(ttl);
    cos_tts = cosd(tts);
    sin_tts = sind(tts);
    cos_tto = cosd(tto);
    sin_tto = sind(tto);
    Cs      = cos_ttl*cos_tts;
    Ss      = sin_ttl*sin_tts;
    Co      = cos_ttl*cos_tto;
    So      = sin_ttl*sin_tto;

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
    elseif tto<90
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

    t1 = 2Cs * Co + Ss * So * cos_psi;
    t2 = 0;
    if bt2 > 0
        t2 = sin(bt2) * ( 2ds * doo + Ss * So * cos(bt1) * cos(bt3) );
    end

    denom = 2 * FT(pi)^2;
    frho  = ((pi-bt2) * t1 + t2) / denom;
    ftau  = (-bt2     * t1 + t2) / denom;

    # fill the values to container
    container[1] = chi_s;
    container[2] = abs(chi_o);
    container[3] = frho;
    container[4] = ftau;

    return nothing
end
