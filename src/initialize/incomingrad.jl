###############################################################################
#
# Create IncomingRadiation
#
###############################################################################
"""
    create_incoming_radiation(FType, swl:Array, file)

Create an `AbstractIncomingRadiation` struct, given
- `FType` Floating number type
- `swl` Standard wave length
- `file` Input file name
"""
function create_incoming_radiation(swl::Array{FT,1}, file::String = file_Sun) where {FT<:AbstractFloat}
    N = length(swl)-1

    # Read data
    _suni  = matread(file)["sun"]
    _wl    = _suni["wl"      ]
    _Edir  = _suni["Edirect" ]
    _Ediff = _suni["Ediffuse"]

    # create arrays
    wl    = zeros(FT, N)
    Edir  = zeros(FT, N)
    Ediff = zeros(FT, N)

    # fill in the arrays
    # println("Reading Optical Parameters from ", swl[1], " to ", swl[end], " length: ", length(swl))
    for i in 1:N
        wo = findall( (_wl.>=swl[i]) .& (_wl.<swl[i+1]) )
        if length(wo)==0
            println("Warning, some wavelengths out of bounds ", swl[i])
        end
        wl[i]    = mean(   _wl[wo])
        Edir[i]  = mean( _Edir[wo])
        Ediff[i] = mean(_Ediff[wo])
    end

    # create struct from the arrays
    return IncomingRadiation{FT}(wl, Edir, Ediff)
end
