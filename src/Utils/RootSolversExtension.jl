module RootSolversExtension

using DocStringExtensions
using RootSolvers

RSM = RootSolvers.RootSolvingMethod
TOL = RootSolvers.AbstractTolerance

export BisectionMethod,
       NewtonBisectionMethod

export find_zero_ext




###############################################################################
#
# New root finding method
#
###############################################################################
"""
    struct BisectionMethod

# Fields
$(FIELDS)
"""
struct BisectionMethod{FT} <: RSM{FT}
    "lower bound"
    x_min::FT
    "upper bound"
    x_max::FT
end




"""
    struct NewtonBisectionMethod

# Fields
$(FIELDS)
"""
struct NewtonBisectionMethod{FT} <: RSM{FT}
    "lower bound"
    x_min::FT
    "upper bound"
    x_max::FT
end








###############################################################################
#
# Judge if root is found
#
###############################################################################
"""
    if_break(tol::TOL{FT}, x1::FT, x2::FT, y::FT)

Determine whether to break, given
- `tol` Tolerance type
- `x1` Lower bound of x
- `x2` Upper bound of x
- `y` Residual of y
"""
function if_break(
            tol::SolutionTolerance{FT},
            x1::FT,
            x2::FT,
            y::FT
            ) where {FT<:AbstractFloat}
    return abs(x1-x2) < tol.tol
end

function if_break(
            tol::ResidualTolerance{FT},
            x1::FT,
            x2::FT,
            y::FT
            ) where {FT<:AbstractFloat}
    return abs(y) < tol.tol
end








###############################################################################
#
# Find the first 0 (root)
#
###############################################################################
"""
    find_zero_ext(f::F, ms::BisectionMethod{FT}, tol::TOL{FT})

Find the solution, given
- `f` A function to solve
- `ms` [`BisectionMethod`](@ref) or [`NewtonBisectionMethod`](@ref)
- `tol` Solution tolerance
"""
function find_zero_ext(
            f::F,
            ms::BisectionMethod{FT},
            tol::TOL{FT}
            ) where {F<:Function, FT<:AbstractFloat}
    # calculate the values for x_min and x_max first
    x_min = ms.x_min;
    x_max = ms.x_max;
    x_mid = x_min   ;
    y_min = f(x_min);
    y_max = f(x_max);

    # if y_min and y_max are on the same side
    if y_min * y_max > 0
        if abs(y_max) < abs(y_min)
            return x_max
        else
            return x_min
        end

    # if y_min is 0
    elseif y_min==0
        return x_min

    # if y_min and y_max are on different side, or one of them is 0
    else
        while true
            x_mid = (x_min + x_max) / 2;
            y_mid = f(x_mid);

            # if difference is lower than the tolerance
            if if_break(tol, x_min, x_max, y_mid)
                break
            end

            # if y_mid is on the same side with one of them
            if y_mid * y_min > 0
                x_min = x_mid;
                y_min = y_mid;
            elseif y_mid * y_max > 0
                x_max = x_mid;
                y_max = y_mid;

            # if y_mid = 0 and one of y_min and y_max is 0
            elseif (y_min==0) && (y_mid==0)
                x_min = x_mid;
                y_min = y_mid;
            elseif (y_max==0) && (y_mid==0)
                x_max = x_mid;
                y_max = y_mid;

            # else set x_max to x_mid (takes the fisrt one)
            else
                x_max = x_mid;
                y_max = y_mid;
            end
        end

        return x_mid
    end
end

function find_zero_ext(
            f::F,
            ms::NewtonBisectionMethod{FT},
            tol::TOL{FT}
            ) where {F<:Function, FT<:AbstractFloat}
    # calculate the values for x_min and x_max first
    x_dif = tol.tol ;
    x_min = ms.x_min;
    x_max = ms.x_max;
    y_min = f(x_min);
    y_max = f(x_max);

    # if y_min and y_max are on the same side
    if y_min * y_max > 0
        if abs(y_max) < abs(y_min)
            return x_max
        else
            return x_min
        end

    # if y_min is 0
    elseif y_min==0
        return x_min

    # if y_min and y_max are on different side, or one of them is 0
    else
        x_ntr = (x_max + x_min) / 2;
        while true
            y_ntr = f(x_ntr);

            # if difference is lower than the tolerance
            if if_break(tol, x_min, x_max, y_ntr)
                break
            end

            # if y_ntr is on the same side with one of them
            if y_ntr * y_min > 0
                x_min = x_ntr;
                y_min = y_ntr;
            elseif y_ntr * y_max > 0
                x_max = x_ntr;
                y_max = y_ntr;

            # if y_ntr = 0 and one of y_min and y_max is 0
            elseif (y_min==0) && (y_ntr==0)
                x_min = x_ntr;
                y_min = y_ntr;
            elseif (y_max==0) && (y_ntr==0)
                x_max = x_ntr;
                y_max = y_ntr;

            # else set x_max to x_ntr (takes the fisrt one)
            else
                x_max = x_ntr;
                y_max = y_ntr;
            end

            #= used for debugging
            @show x_ntr;
            @show y_ntr;
            @show x_min;
            @show x_max;
            @show x_dif;
            println("");
            sleep(0.1)
            =#

            # update x_ntr using Newton Raphson
            x_dx   = x_ntr + x_dif * 10;
            y_dx   = f(x_dx);
            slope  = (y_dx - y_ntr) / (x_dx - x_ntr);
            x_ntr -= y_ntr / slope;

            # set x_ntr in the correct range
            if (x_ntr>=x_max) || (x_ntr<=x_min) || isnan(x_ntr)
                x_ntr = (x_min + x_max) / 2;
            end
        end

        return x_ntr
    end
end




end # module
