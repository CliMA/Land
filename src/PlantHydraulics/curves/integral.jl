function vc_integral(
            vc::LogisticSingle{FT},
            p_dos::FT,
            p_ups::FT
) where {FT<:AbstractFloat}
    @assert p_ups <= 0 && p_dos <= 0;

    # unpack data from VC
    @unpack a,b = vc;
    _t_dos = log(a * exp(b*p_dos) + 1);
    _t_ups = log(a * exp(b*p_ups) + 1);

    return (a+1) / (a*b) * (_t_ups - _t_dos)
end




function vc_integral(
            vc::LogisticSingle{FT},
            p_dos::FT,
            p_ups::FT,
            h::FT,
            E::FT,
            Kmax::FT
) where {FT<:AbstractFloat}
    @assert p_ups <= 0 && p_dos <= 0;

    # unpack data from VC
    @unpack a,b = vc;
    _krghe = Kmax * ρg_MPa(FT) * h * (a+1) / a + E;
    _lower = b * _krghe;
    _multi = Kmax * (a+1) / a * E;
    _upper_dos = log(a * _krghe * exp(b*p_dos) + E);
    _upper_ups = log(a * _krghe * exp(b*p_ups) + E);

    return _multi * (_upper_ups - _upper_dos) / _lower
end




function vc_integral(
            vc::PowerSingle{FT},
            p_dos::FT,
            p_ups::FT
) where {FT<:AbstractFloat}
    @assert p_ups <= 0 && p_dos <= 0;

    # unpack data from vc
    @unpack a,b = vc;
    _f_dos = p_dos * _₂F₁(1, 1/b, 1+1/b, -a*(-p_dos)^b);
    _f_ups = p_ups * _₂F₁(1, 1/b, 1+1/b, -a*(-p_ups)^b);

    return _f_ups - _f_dos
end




function vc_integral(
            vc::PowerSingle{FT},
            p_dos::FT,
            p_ups::FT,
            h::FT,
            E::FT,
            Kmax::FT
) where {FT<:AbstractFloat}
    @assert p_ups <= 0 && p_dos <= 0;

    # unpack data from VC
    @unpack a,b = vc;
    _krghe = Kmax * ρg_MPa(FT) * h + E;
    _multi = Kmax * E / _krghe;
    _f_dos = p_dos * _₂F₁(1, 1/b, 1+1/b, -E * a * (-p_dos)^b / _krghe);
    _f_ups = p_ups * _₂F₁(1, 1/b, 1+1/b, -E * a * (-p_ups)^b / _krghe);

    return _multi * (_f_ups - _f_dos)
end




function vc_integral(
            vc::WeibullSingle{FT},
            p_dos::FT,
            p_ups::FT
) where {FT<:AbstractFloat}
    @assert p_ups <= 0 && p_dos <= 0;

    # unpack data from VC
    @unpack b,c = vc;

    # compute the incomplete gamma function
    _Γ_dos = gamma(1/c, (-p_dos/b)^c);
    _Γ_ups = gamma(1/c, (-p_ups/b)^c);

    return b/c * (_Γ_ups - _Γ_dos)
end
