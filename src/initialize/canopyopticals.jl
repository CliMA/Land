###############################################################################
#
# Create CanopyOpticals
#
###############################################################################
"""
create_canopy_opticals(FT, nWL::Int, nLayer::Int, nAzi::Int, nIncl::Int)

Create a canopy optical properties struct, given
- `FT` Floating number type
- `nWL` Number of wave length
- `nLayer` Number of canopy layers
- `nAzi` Number of arimuth angles
- `AIncl` Number of inclination angles
"""
function create_canopy_opticals(
            FT,
            rt_dim::RTDimentions
)
    @unpack nAzi, nIncl, nLayer, nWL = rt_dim;

    return CanopyOpticals{FT}(nAzi   = nAzi  ,
                              nIncl  = nIncl ,
                              nLayer = nLayer,
                              nWL    = nWL   )
end
