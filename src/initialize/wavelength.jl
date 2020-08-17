###############################################################################
#
# Initialize WaveLengths
#
###############################################################################
"""
    create_wave_length(FT)

Create [`WaveLengths`](@ref) type struct, given
- `FT` Floating number type
"""
function create_wave_length(FT)
    sWL::Array{FT,1} = [collect(FT(400.0):FT(10.0):FT( 650.1));
                        collect(FT(655.0):FT( 5.0):FT( 770.1));
                        collect(FT(780.0):FT(25.0):FT(2400.1))]
    optis = create_leaf_opticals(sWL, FILE_OPTI);

    return WaveLengths{FT}(sWL=sWL, optis=optis)
end
