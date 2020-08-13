
function create_canopy_rt(
            FT;
            nLayer::Int = 20
)
    return Canopy4RT{FT}(nLayer = nLayer)
end
