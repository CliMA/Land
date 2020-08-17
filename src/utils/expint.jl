###############################################################################
#
# ?
#
###############################################################################
"""
    expint(x::FT) where {FT<:AbstractFloat}

# TODO add function description
"""
function expint(
            x::FT
) where {FT<:AbstractFloat}
    pn = Polynomial(FT[ 8.267661952366478e+00,
                       -7.773807325735529e-01,
                       -3.012432892762715e-01,
                       -7.811863559248197e-02,
                       -1.019573529845792e-02,
                       -6.973790859534190e-04,
                       -2.569498322115933e-05,
                       -4.819538452140960e-07,
                       -3.602693626336023e-09]);
    egamma = FT(0.57721566490153286061);
    polyv  = pn(real(x));
    if abs(imag(x)) <= polyv
        #initialization
        xk = x;
        yk = -egamma - log(xk);
        j = 1;
        pterm = xk;
        term = xk;
        while abs(term) > (eps(yk))
            yk = yk + term;
            j = j + 1;
            pterm = -xk.*pterm/j;
            term = pterm/j;
        end # end of the while loop
        y = yk;
    else
        n = FT(1);
        xk = x;
        am2 = FT(0)
        bm2 = FT(1)
        am1 = FT(1)
        bm1 = xk
        f = am1 / bm1;
        oldf = Inf;
        j = FT(2);

        while abs(f-oldf) > (100*eps(FT)*abs(f))
            alpha = n-1+(j/2); # note: beta= 1
            #calculate A(j), B(j), and f(j)
            a = am1 + alpha * am2;
            b = bm1 + alpha * bm2;

            # save new normalized variables for next pass through the loop
            #  note: normalization to avoid overflow or underflow
            am2 = am1 / b;
            bm2 = bm1 / b;
            am1 = a / b;
            bm1 = FT(1);

            f = am1;
            j = j+1;

            # calculate the coefficients for j odd
            alpha = (j-1)/2;
            beta = xk;
            a = beta * am1 + alpha * am2;
            b = beta * bm1 + alpha * bm2;
            am2 = am1 / b;
            bm2 = bm1 / b;
            am1 = a / b;
            bm1 = 1;
            oldf = f;
            f = am1;
            #println(typeof(xk), typeof(f), typeof(alpha))
            j = j+1;
        end

        y= exp(-xk) * f - 1im*FT(pi)*((real(xk)<0)&(imag(xk)==0));
    end

    return y
end
