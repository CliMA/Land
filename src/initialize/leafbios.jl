###############################################################################
#
# Create LeafBios
#
###############################################################################
"""
    create_leaf_bios(FT, rt_dim::RTDimensions)

Create a [`LeafBios`](@ref) type struct, given
- `FT` Floating number type
- `rt_dim` [`RTDimensions`](@ref) type struct
"""
function create_leaf_bios(FT, rt_dim::RTDimensions)
    @unpack nWL, nWLE, nWLF = rt_dim;

    return LeafBios{FT}(nWL  = nWL ,
                        nWLE = nWLE,
                        nWLF = nWLF)
end
