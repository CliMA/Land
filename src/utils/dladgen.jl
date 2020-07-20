#=
"""
    dcum(a::FT, b::FT, t::FT)

TODO Add function description
"""
=#
function dcum(a::FT, b::FT, t::FT) where {FT}
    y = FT(0.0)
    if a >= 1
        f = 1 - cosd(t)
    else
        epsi = 1e-8
        delx = FT(1.0)
        x    = 2 * deg2rad(t)
        p    = x
        # add a iteration number (n) control to avoid rounding issue in Float32
        n    = 0
        while (delx >= epsi) && (n<50)
            n   += 1
            y    = a*sin(x) + b*sin(2x)/2
            dx   = (y - x + p) / 2
            x   += dx
            delx = abs(dx)
        end
    	f = (2*y + p) / pi
    end
    return f
end




#=
"""
    dladgen(a::FT, b::FT, litab_bnd::Array)

TODO Calculate the freqency of WHAT?
"""
=#
function dladgen(a::FT, b::FT, litab_bnd::Array{FT,2}) where {FT}
    @assert a+b< 1
    freq = similar(litab_bnd[:,1])
    for i in eachindex(freq)
        freq[i] = dcum(a, b, litab_bnd[i,2]) - dcum(a, b, litab_bnd[i,1])
    end
    return freq
end
