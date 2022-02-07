#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2021-Aug-10: refactor the structure with renamed fields
#     2021-Aug-10: add a constructor within the structure to avoid external initialization
#     2021-Oct-19: sort variable to prognostic and dignostic catergories
#     2022-Jan-24: fix documentation
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores wave length information.

# Fields

$(TYPEDFIELDS)

"""
mutable struct WaveLengthSet{FT<:AbstractFloat}
    # parameters that do not change with time
    "Indicies of Λ_NIR in Λ"
    IΛ_NIR::Vector{Int}
    "Indicies of Λ_PAR in Λ"
    IΛ_PAR::Vector{Int}
    "Indicies of Λ_SIF in Λ"
    IΛ_SIF::Vector{Int}
    "Indicies of Λ_SIFE in Λ"
    IΛ_SIFE::Vector{Int}
    "Number of wavelength bins"
    NΛ::Int
    "Number of wavelength bins for PAR"
    NΛ_PAR::Int
    "Number of wavelength bins for SIF"
    NΛ_SIF::Int
    "Number of wavelength bins for SIF excitation"
    NΛ_SIFE::Int
    "Wavelength limits for NIR `[nm]`"
    WL_NIR::Vector{FT}
    "Wavelength limits for PAR `[nm]`"
    WL_PAR::Vector{FT}
    "Wavelength limits for SIF emission `[nm]`"
    WL_SIF::Vector{FT}
    "Wavelength limits for SIF excitation `[nm]`"
    WL_SIFE::Vector{FT}
    "Standard wavelength (boundaries) `[nm]`"
    SΛ::Vector{FT}
    "Differential wavelength `[nm]`"
    ΔΛ::Vector{FT}
    "Differential wavelength for PAR `[nm]`"
    ΔΛ_PAR::Vector{FT}
    "Differential wavelength for SIF excitation `[nm]`"
    ΔΛ_SIFE::Vector{FT}
    "Wavelength (bins) `[nm]`"
    Λ::Vector{FT}
    "Wavelength bins for PAR `[nm]`"
    Λ_PAR::Vector{FT}
    "Wavelength bins for SIF `[nm]`"
    Λ_SIF::Vector{FT}
    "Wavelength bins for SIF excitation `[nm]`"
    Λ_SIFE::Vector{FT}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2021-Nov-24: isolate the constructor from structure
#     2022-Jan-24: fix documentation
#
#######################################################################################################################################################################################################
"""

    WaveLengthSet{FT}(swl::Vector=WAVELENGTHS; opti::String=OPTI_2021) where {FT<:AbstractFloat}

Constructor for [`WaveLengthSet`](@ref), given
- `swl` Standard wave length boundaries `[nm]`, default is `ClimaCache.WAVELENGTHS`
- `opti` Optical properties file path, default is `ClimaCache.OPTI_2021`

---
# Examples
```julia
wls = WaveLengthSet{Float64}();
wls = WaveLengthSet{Float64}(collect(400:5:2500));
wls = WaveLengthSet{Float64}(collect(400:5:2500); opti=ClimaCache.OPTI_2017);
```
"""
WaveLengthSet{FT}(swl::Vector=WAVELENGTHS; opti::String=OPTI_2021) where {FT<:AbstractFloat} = (
    _dwl    = diff(swl);
    _λ      = zeros(FT, length(swl)-1);
    _opti   = matread(opti)["optipar"];
    _λ_opti = _opti["wl"];
    @inbounds for _i in 1:length(swl)-1
        _wo    = findall( swl[_i] .<= _λ_opti .< swl[_i+1] );
        _λ[_i] = mean(_λ_opti[_wo]);
    end;
    _iλ_nir  = findall( 700 .<= _λ .<= 2500 );
    _iλ_par  = findall( 400 .<= _λ .<= 700  );
    _iλ_sif  = findall( 640 .<= _λ .<= 850  );
    _iλ_sife = findall( 400 .<= _λ .<= 750  );

    return WaveLengthSet{FT}(_iλ_nir, _iλ_par, _iλ_sif, _iλ_sife, length(_λ), length(_iλ_par), length(_iλ_sif), length(_iλ_sife), FT[700,2500], FT[400,700], FT[640,850], FT[400,750], swl, _dwl,
                             _dwl[_iλ_par], _dwl[_iλ_sife], _λ, _λ[_iλ_par], _λ[_iλ_sif], _λ[_iλ_sife])
);
