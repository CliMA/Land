###############################################################################
#
# Initialize WaveLengths
#
###############################################################################
"""
    create_wave_length(
                FT,
                sWLs = [collect(400.0:10.0: 650.1);
                        collect(655.0: 5.0: 770.1);
                        collect(780.0:25.0:2400.1)])

Create [`WaveLengths`](@ref) type struct, given
- `FT` Floating number type
"""
function create_wave_length(
            FT,
            sWLs = [collect(400.0:10.0: 650.1);
                    collect(655.0: 5.0: 770.1);
                    collect(780.0:25.0:2400.1)]
)
    sWL   = FT.(sWLs);
    optis = create_leaf_opticals(sWL, OPTI_2021);

    return WaveLengths{FT}(sWL=sWL, optis=optis)
end
