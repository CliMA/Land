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
function create_canopy_opticals(FT, nWL::Int, nLayer::Int, nAzi::Int, nIncl::Int)
    return CanopyOpticals{FT}(nWL=nWL, nLayer=nLayer, nAzi=nAzi, nIncl=nIncl)
end
