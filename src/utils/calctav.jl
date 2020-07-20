###############################################################################
#
# Compute transmission of isotropic radiation
#
###############################################################################
"""
    calctav(α::FT, nr::FT) where {FT<:AbstractFloat}

Computes transmission of isotropic radiation across an interface between two dielectrics (Stern F. (1964), Allen W.A. (1973)). From calctav.m in PROSPECT-D
- `α` angle of incidence
- `nr` Index of refraction
"""
function calctav(
            α::FT,
            nr::FT
) where {FT<:AbstractFloat}
    tt = typeof(nr)
    rd  = tt(pi/180)
    n2  = nr^2
    np  = n2+1
    nm  = n2-1
    a   = (nr+1)*(nr+1)/2
    k   = -(n2-1)*(n2-1)/4
    sa  = sin(α*rd)

    if α!=90
        b1  = sqrt((sa^2-np/2).*(sa^2-np/2)+k)
    else
        b1 = 0
    end

    b2  = sa^2-np/2
    b   = b1-b2
    b3  = b^3
    a3  = a^3
    ts  = (k^2/(6*b3)+k/b-b/2)-(k^2/(6*a3)+k/a-a/2)
    tp1 = -2*n2*(b-a)/(np^2)
    tp2 = -2*n2*np*log(b/a)/(nm^2)
    tp3 = n2*(1/b-1/a)/2
    tp4 = 16*n2^2*(n2^2+1)*log((2*np*b-nm^2)/(2*np*a-nm^2))/(np^3*nm^2)
    tp5 = 16*n2^3*(1/(2*np*b-nm^2)-1/(2*np*a-nm^2))/(np^3)
    tp  = tp1+tp2+tp3+tp4+tp5
    tav = (ts+tp)/(2*sa^2)

    return tav
end
