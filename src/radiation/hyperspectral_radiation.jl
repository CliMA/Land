#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2021-Oct-22: refactor the structure with renamed fields
#     2021-Oct-22: add a constructor to define the structure from wavelength sets and prescribed wave shape
#
#######################################################################################################################################################################################################
"""
$(TYPEDEF)

Structure that stores hyperspectral radiation information

# Fields
$(TYPEDFIELDS)
"""
mutable struct HyperspectralRadiation{FT<:AbstractFloat}
    # prognostic variables that change with time
    "Direct radiation `[mW m⁻² nm⁻¹]`"
    e_direct::Vector{FT}
    "Diffuse radiation `[mW m⁻² nm⁻¹]`"
    e_diffuse::Vector{FT}
end


"""
    HyperspectralRadiation{FT}(wls::WaveLengthSet = WaveLengthSet{FT}(), wlfn::String = FILE_SUN)

Constructor for [`HyperspectralRadiation`](@ref), given
- `wls` [`WaveLengthSet`](@ref) type struct that defines wavelength settings
- `wlfn` File path to solar radiation setting, default is `ClimaCache.FILE_SUN`

---
# Examples
```julia
rad = HyperspectralRadiation{FT}();
rad = HyperspectralRadiation{FT}(WaveLengthSet{FT}(collect(400:50:2400)));
rad = HyperspectralRadiation{FT}(WaveLengthSet{FT}(collect(400:50:2400)), "");
rad = HyperspectralRadiation{FT}(WaveLengthSet{FT}(collect(400:50:2400)), ClimaCache.FILE_SUN);
```
"""
HyperspectralRadiation{FT}(wls::WaveLengthSet = WaveLengthSet{FT}(), wlfn::String = FILE_SUN) where {FT<:AbstractFloat} = (
    @unpack SΛ, NΛ = wls;

    # create arrays
    _e_direct  = zeros(FT, NΛ);
    _e_diffuse = zeros(FT, NΛ);

    # Read data from SCOPE MAT file if file is given, if not use 0
    if isfile(wlfn)
        _suni    = matread(wlfn)["sun"];
        _wl      = _suni["wl"      ];
        __e_dir  = _suni["Edirect" ];
        __e_diff = _suni["Ediffuse"];

        # fill in the arrays
        for _i in 1:NΛ
            _wi = findall( (_wl.>=SΛ[_i]) .& (_wl.<SΛ[_i+1]) );
            if length(_wi)==0 @warn "Some wavelengths out of bounds $(string(SΛ[_i]))" end;
            _e_direct[_i]  = mean( __e_dir[_wi]);
            _e_diffuse[_i] = mean(__e_diff[_wi]);
        end;
    end;

    return HyperspectralRadiation{FT}(_e_direct, _e_diffuse)
);
