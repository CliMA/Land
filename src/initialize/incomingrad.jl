###############################################################################
#
# Create IncomingRadiation
#
###############################################################################
"""
    create_incoming_radiation(sWL:Array, file)

Create an `AbstractIncomingRadiation` struct, given
- `FType` Floating number type
- `sWL` Standard wave length
- `file` Input file name
"""
function create_incoming_radiation(
            wl_set::WaveLengths{FT},
            file::String = file_Sun
) where {FT<:AbstractFloat}
    @unpack sWL, nWL = wl_set;

    # Read data
    _suni  = matread(file)["sun"]
    _wl    = _suni["wl"      ]
    _Edir  = _suni["Edirect" ]
    _Ediff = _suni["Ediffuse"]

    # create arrays
    # WL  = zeros(FT, N)
    Edir  = zeros(FT, nWL)
    Ediff = zeros(FT, nWL)

    # fill in the arrays
    # println("Reading Optical Parameters from ", sWL[1], " to ", sWL[end], " length: ", length(sWL))
    for i in 1:nWL
        wo = findall( (_wl.>=sWL[i]) .& (_wl.<sWL[i+1]) )
        if length(wo)==0
            println("Warning, some wavelengths out of bounds ", sWL[i])
        end
        # remove WL to avoid duplication
        # WL[i]  = mean(   _wl[wo])
        Edir[i]  = mean( _Edir[wo])
        Ediff[i] = mean(_Ediff[wo])
    end

    # create struct from the arrays
    return IncomingRadiation{FT}(Edir, Ediff)
end
