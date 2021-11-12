###############################################################################
#
# Xylem K ratio functions
#
###############################################################################
"""
    xylem_k_ratio(
                vc::AbstractXylemVC{FT},
                p_25::FT,
                vis::FT
    ) where {FT<:AbstractFloat}

Returns the relative hydraulic conductance, given
- `vc` Xylem vulnerability curve
- `p` Xylem pressure at 298.15 K in `[MPa]`
- `p_25` Equivalent xylem pressure at 298.15 K in `[MPa]`
- `vis` Relative viscosity. If missing, vis = 1.
"""
function xylem_k_ratio(
            vc::LogisticSingle{FT},
            p_25::FT,
            vis::FT = FT(1)
) where {FT<:AbstractFloat}
    if p_25>=0
        return 1 / vis
    end

    @unpack a,b = vc;
    return max( FT(1e-4), (1 - 1/(1 + a * exp(b * p_25))) * (a+1)/a ) / vis
end




function xylem_k_ratio(
            vc::PowerSingle{FT},
            p_25::FT,
            vis::FT = FT(1)
) where {FT<:AbstractFloat}
    if p_25>=0
        return 1 / vis
    end

    @unpack a,b = vc;
    return max( FT(1e-4), 1 / (1 + a*(-p_25)^b) ) / vis
end




function xylem_k_ratio(
            vc::WeibullSingle{FT},
            p_25::FT,
            vis::FT = FT(1)
) where {FT<:AbstractFloat}
    @unpack b,c = vc;

    if p_25<0
        kr = max( FT(1e-4), exp( -1 * (-p_25 / b) ^ c ) ) / vis;
    else
        kr = 1 / vis;
    end

    return kr
end




function xylem_k_ratio(
            vc::WeibullDual{FT},
            p_25::FT,
            vis::FT = FT(1)
) where {FT<:AbstractFloat}
    @unpack b1,c1,f1,b2,c2,f2 = vc;

    if p_25<0
        k1 = exp( -1 * (-p_25 / b1) ^ c1 ) * f1;
        k2 = exp( -1 * (-p_25 / b2) ^ c2 ) * f2;
        kr = max( FT(1e-4), (k1+k2) ) / vis;
    else
        kr = 1 / vis;
    end

    return kr
end








###############################################################################
#
# Xylem P_crit functions
#
###############################################################################
"""
    xylem_p_crit(vc::AbstractXylemVC{FT}, f_st::FT) where {FT<:AbstractFloat}

Returns the relative hydraulic conductance, given
- `vc` Xylem vulnerability curve
- `st` Relative surface tension. If missing, vis = 1.
"""
function xylem_p_crit(
            vc::WeibullSingle{FT},
            f_st::FT = FT(1)
) where {FT<:AbstractFloat}
    @unpack b,c = vc;

    return -b * log( FT(1000) ) ^ (1 / c) * f_st
end




function xylem_p_crit(
            vc::WeibullDual{FT},
            f_st::FT = FT(1)
) where {FT<:AbstractFloat}
    @unpack b1,c1,b2,c2 = vc;

    p1 = -b1 * log( FT(1000) ) ^ (1 / c1) * f_st
    p2 = -b2 * log( FT(1000) ) ^ (1 / c2) * f_st

    return min(p1, p2)
end








###############################################################################
#
# Fit xylem vulnerability curve
#
###############################################################################

function fit_xylem_VC(xs::Array{FT,1}, ys::Array{FT,1}; label="TPLC") where {FT<:AbstractFloat}
    _ps = xs;
    _ks = ys;

    @inline f(x) = (
        @show x;
        _vc = WeibullSingle{FT}(b=x[1], c=x[2]);
        _kp = xylem_k_ratio.([_vc], _ps);
        return -sum((_kp .- _ks) .^ 2);
    );

    @inline g(x) = (
        @show x;
        _vc = WeibullSingle{FT}(b=x[1], c=x[2]);
        _kp = x[3] * xylem_k_ratio.([_vc], _ps);
        return -sum((_kp .- _ks) .^ 2);
    );

    # if x is tension and y is plc
    if label=="TPLC"
        _ps = -xs;
        _ks = (100 .- ys) ./ 100;

        _st = SolutionToleranceND{FT}([1e-3, 1e-3], 30);
        _ms = ReduceStepMethodND{FT}(x_mins = FT[1e-3, 1e-3],
                                     x_maxs = FT[ 100, 100],
                                     x_inis = [1, 1],
                                     Δ_inis = FT[1, 1]);
        _bc = find_peak(f, _ms, _st);

        @show _bc;

        return _bc
    elseif label=="PK"
        _ps = xs;
        _ks = ys;

        _st = SolutionToleranceND{FT}([1e-3, 1e-3, 1e-2], 30);
        _ms = ReduceStepMethodND{FT}(x_mins = FT[1e-3, 1e-3, 1e-3],
                                     x_maxs = FT[ 100, 100, 100],
                                     x_inis = [1, 1, 50],
                                     Δ_inis = FT[1, 1, 10]);
        _bck = find_peak(g, _ms, _st);

        @show _bck;

        return _bck
    else
        return nothing
    end
end
