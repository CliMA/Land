#=
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
        _vc = WeibullVC{FT}(x[1], x[2]);
        _kp = relative_hydraulic_conductance.([_vc], _ps);
        return -sum((_kp .- _ks) .^ 2);
    );

    @inline g(x) = (
        @show x;
        _vc = WeibullVC{FT}(x[1], x[2]);
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
=#
