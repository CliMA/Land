###############################################################################
#
# Create IncomingRadiation
#
###############################################################################
"""
    create_incoming_radiation(
                wls::WaveLengths{FT},
                wlfn::String = FILE_SUN
    ) where {FT<:AbstractFloat}

Create an `AbstractIncomingRadiation` struct, given
- `wls` [`WaveLengths`](@ref) type struct
- `wlfn` File that saves incoming wave information
"""
function create_incoming_radiation(
            wls::WaveLengths{FT},
            wlfn::String = FILE_SUN
) where {FT<:AbstractFloat}
    @unpack sWL, nWL = wls;

    # Read data
    _suni  = matread(wlfn)["sun"]
    _wl    = _suni["wl"      ]
    _Edir  = _suni["Edirect" ]
    _Ediff = _suni["Ediffuse"]

    # create arrays
    # WL  = zeros(FT, N)
    Edir  = zeros(FT, nWL)
    Ediff = zeros(FT, nWL)

    # fill in the arrays
    # println("Reading Optical Parameters from ", sWL[1], " to ", sWL[end])
    for i in 1:nWL
        wo = findall( (_wl.>=sWL[i]) .& (_wl.<sWL[i+1]) )
        if length(wo)==0
            @warn "Some wavelengths out of bounds $(string(sWL[i]))";
        end
        # remove WL to avoid duplication
        # WL[i]  = mean(   _wl[wo])
        Edir[i]  = mean( _Edir[wo])
        Ediff[i] = mean(_Ediff[wo])
    end

    # create struct from the arrays
    return IncomingRadiation{FT}(Edir, Ediff)
end
