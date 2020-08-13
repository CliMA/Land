###############################################################################
#
# Create LeafBios
#
###############################################################################
"""
    create_leaf_bios(FT, nWL::Int, nWLE::Int, nWLF::Int)

Create a leaf biological parameters struct, given
- `FType` Floating number type
- `nWL` Number of wave length
- `nWLE` Number of excitation wave length
- `nWLF` Number of fluorescence wave length

Returns a [`LeafBios`](@ref) type struct.
"""
function create_leaf_bios(FT, nWL::Int, nWLE::Int, nWLF::Int)
    return LeafBios{FT}(nWL=nWL, nWLE=nWLE, nWLF=nWLF)
end
