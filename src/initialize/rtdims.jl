###############################################################################
#
# Initialize RTDimensions
#
###############################################################################
"""
    create_rt_dims(
                can::Canopy4RT{FT},
                wls::WaveLengths{FT}
    ) where {FT<:AbstractFloat}

Create [`RTDimensions`](@ref), given
- `can` [`Canopy4RT`](@ref) type struct
- `wls` [`WaveLengths`](@ref) type struct
"""
function create_rt_dims(
            can::Canopy4RT{FT},
            wls::WaveLengths{FT}
) where {FT<:AbstractFloat}
    @unpack nAzi, nIncl, nLayer = can;
    @unpack nPAR, nWL, nWLE, nWLF = wls;

    return RTDimensions(nAzi   = nAzi  ,
                        nIncl  = nIncl ,
                        nLayer = nLayer,
                        nPAR   = nPAR  ,
                        nWL    = nWL   ,
                        nWLE   = nWLE  ,
                        nWLF   = nWLF  )
end
