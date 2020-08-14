
function create_rt_dims(
            can_rt::Canopy4RT{FT},
            wl_set::WaveLengths{FT}
) where {FT<:AbstractFloat}
    @unpack nAzi, nIncl, nLayer = can_rt;
    @unpack nPAR, nWL, nWLE, nWLF = wl_set;

    return RTDimentions(nAzi   = nAzi  ,
                        nIncl  = nIncl ,
                        nLayer = nLayer,
                        nPAR   = nPAR  ,
                        nWL    = nWL   ,
                        nWLE   = nWLE  ,
                        nWLF   = nWLF  )
end
