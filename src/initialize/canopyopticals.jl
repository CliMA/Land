###############################################################################
#
# Initialize CanopyOpticals
#
###############################################################################
"""
    create_canopy_opticals( FT, rt_dim::RTDimensions)

Create a [`CanopyOpticals`](@ref) struct, given
- `FT` Floating number type
- `rt_dim` [`RTDimensions`](@ref) type struct
"""
function create_canopy_opticals(FT, rt_dim::RTDimensions)
    @unpack nAzi, nIncl, nLayer, nWL = rt_dim;

    return CanopyOpticals{FT}(nAzi   = nAzi  ,
                              nIncl  = nIncl ,
                              nLayer = nLayer,
                              nWL    = nWL   )
end
