###############################################################################
#
# Initialize CanopyOpticals
#
###############################################################################
"""
    create_canopy_rads(FT, rt_dim::RTDimensions)

Create a [`CanopyRads`](@ref) struct, given
- `FT` Floating number type
- `rt_dim` [`RTDimensions`](@ref) type struct
"""
function create_canopy_rads(FT, rt_dim::RTDimensions)
    @unpack nAzi, nIncl, nLayer, nLevel, nWL, nWLF = rt_dim;

    return CanopyRads{FT}(nAzi   = nAzi  ,
                          nIncl  = nIncl ,
                          nLayer = nLayer,
                          nLevel = nLevel,
                          nWL    = nWL   ,
                          nWLF   = nWLF  )
end
