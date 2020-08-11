###############################################################################
#
# Create LeafBios
#
###############################################################################
"""
    create_leaf_bios(FT, nWl::Int, nWle::Int, nWlf::Int)

Create a leaf biological parameters struct, given
- `FType` Floating number type
- `nWl` Number of wave length
- `nWle` Number of excitation wave length
- `nWlf` Number of fluorescence wave length

Returns a [`LeafBios`](@ref) type struct.
"""
function create_leaf_bios(FT, nWL::Int, nWLe::Int, nWLf::Int)
    return LeafBios{FT}(nWL=nWL, nWLe=nWLe, nWLf=nWLf)
end
