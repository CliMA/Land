###############################################################################
#
# Initialize Canopy4RT
#
###############################################################################
"""
    create_canopy_rt(FT; nLayer::Int = 20, LAI::Number = FT(3))

Create [`Canopy4RT`](@ref), given
- `FT` Floating number type
- `nLayer` Number of canopy layers
- `LAI` Leaf area index
"""
function create_canopy_rt(FT; nLayer::Int = 20, LAI::Number = FT(3))
    return Canopy4RT{FT}(LAI    = LAI   ,
                         nLayer = nLayer)
end
