"""

    create_wave_length(
                FT,
                sWLs = [collect(400.0:10.0: 650.1); collect(655.0: 5.0: 770.1); collect(780.0:25.0:2400.1)];
                max_NIR::Number = 2500,
                max_PAR::Number = 700,
                max_SIF::Number = 750,
                min_NIR::Number = 700,
                min_PAR::Number = 400,
                min_SIF::Number = 400,
                opti_file::String = OPTI_2021)

Create [`WaveLengths`](@ref) type struct, given
- `FT` Floating number type
- `sWLs` Shortwave wavelength bins
- `max_NIR` Maximal NIR wavelength
- `max_PAR` Maximal PAR wavelength
- `max_ISF` Maximal SIF excitation wavelength
- `min_NIR` Minimal NIR wavelength
- `min_PAR` Minimal PAR wavelength
- `min_ISF` Minimal SIF excitation wavelength
- `opti_file` Input reference optical file path

"""
function create_wave_length(
            FT,
            sWLs = [collect(400.0:10.0: 650.1); collect(655.0: 5.0: 770.1); collect(780.0:25.0:2400.1)];
            max_NIR::Number = 2500,
            max_PAR::Number = 700,
            max_SIF::Number = 750,
            min_NIR::Number = 700,
            min_PAR::Number = 400,
            min_SIF::Number = 400,
            opti_file::String = OPTI_2021)
    sWL   = FT.(sWLs);
    optis = create_leaf_opticals(sWL, opti_file);

    return WaveLengths{FT}(minwlPAR = min_PAR, maxwlPAR = max_PAR, minwlNIR = min_NIR, maxwlNIR = max_NIR, minwle = min_SIF, maxwle = max_SIF, sWL = sWL, optis = optis)
end
