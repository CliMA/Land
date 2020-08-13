###############################################################################
#
# Create CanopyOpticals
#
###############################################################################
"""
    create_canopy_rads(FT, nWL::Int, nLayer::Int, nAzi::Int, nIncl::Int)

Create a canopy radiaitons struct, given
- `FT` Floating number type
- `nWL` Number of wave length
- `nLayer` Number of canopy layers
- `nAzi` Number of arimuth angles
- `AIncl` Number of inclination angles
"""
function create_canopy_rads(
            FT,
            rt_dim::RTDimentions
)
    @unpack nAzi, nIncl, nLayer, nLevel, nWL, nWLF = rt_dim;

    return CanopyRads{FT}(nAzi   = nAzi  ,
                          nIncl  = nIncl ,
                          nLayer = nLayer,
                          nLevel = nLevel,
                          nWL    = nWL   ,
                          nWLF   = nWLF  )
end
