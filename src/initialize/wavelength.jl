###############################################################################
#
# Create wavelength
#
###############################################################################
"""
"""
function create_wave_length(FT)
    swl::Array{FT,1} = [collect(FT(400.0):FT(10.0):FT( 650.1));
                        collect(FT(655.0):FT( 5.0):FT( 770.1));
                        collect(FT(780.0):FT(25.0):FT(2400.1))]
    optis = create_leaf_opticals(swl, file_Opti);

    return WaveLengths{FT}(swl=swl, optis=optis)
end
