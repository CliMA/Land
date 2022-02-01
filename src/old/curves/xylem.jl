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
        _kp = relative_hydraulic_conductance.([_vc], _ps);
        return -sum((_kp .- _ks) .^ 2);
    );

    @inline g(x) = (
        @show x;
        _vc = WeibullSingle{FT}(b=x[1], c=x[2]);
        _kp = x[3] * relative_hydraulic_conductance.([_vc], _ps);
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
