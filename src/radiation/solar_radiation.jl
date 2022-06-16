#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jun-15: add abstract type for incoming solar radiation
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractRadiation:
- [`BroadbandRadiation`](@ref)
- [`HyperspectralRadiation`](@ref)
"""
abstract type AbstractRadiation{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-15: add broadband solar radiation
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores broadband radiation information

# Fields

$(TYPEDFIELDS)

"""
mutable struct BroadbandRadiation{FT} <: AbstractRadiation{FT}
    # prognostic variables that change with time
    "Diffuse radiation from NIR region `[W m⁻²]`"
    e_diffuse_nir::FT
    "Diffuse radiation from PAR region `[W m⁻²]`"
    e_diffuse_par::FT
    "Direct radiation from NIR region `[W m⁻²]`"
    e_direct_nir::FT
    "Direct radiation from PAR region `[W m⁻²]`"
    e_direct_par::FT
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2021-Oct-22: add constructor
#
#######################################################################################################################################################################################################
"""

    BroadbandRadiation{FT}() where {FT<:AbstractFloat}

Constructor for [`BroadbandRadiation`](@ref)
"""
BroadbandRadiation{FT}() where {FT<:AbstractFloat} = BroadbandRadiation{FT}(0, 0, 0, 0);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2021-Oct-22: refactor the structure with renamed fields
#     2021-Oct-22: add a constructor to define the structure from wavelength sets and prescribed wave shape
#     2022-Jun-15: change the order of variables
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores hyperspectral radiation information

# Fields

$(TYPEDFIELDS)

"""
mutable struct HyperspectralRadiation{FT} <: AbstractRadiation{FT}
    # prognostic variables that change with time
    "Diffuse radiation `[mW m⁻² nm⁻¹]`"
    e_diffuse::Vector{FT}
    "Direct radiation `[mW m⁻² nm⁻¹]`"
    e_direct::Vector{FT}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2021-Oct-22: add constructor
#     2022-May-25: fix documentation
#     2022-Jun-15: change the order of variables
#
#######################################################################################################################################################################################################
"""

    HyperspectralRadiation{FT}(wls::WaveLengthSet = WaveLengthSet{FT}(); file::String = FILE_SUN) where {FT<:AbstractFloat}

Constructor for [`HyperspectralRadiation`](@ref), given
- `wls` [`WaveLengthSet`](@ref) type struct that defines wavelength settings
- `file` File path to solar radiation setting, default is `ClimaCache.FILE_SUN`

---
# Examples
```julia
rad = HyperspectralRadiation{FT}();
rad = HyperspectralRadiation{FT}(WaveLengthSet{FT}(collect(400:50:2400)));
rad = HyperspectralRadiation{FT}(WaveLengthSet{FT}(collect(400:50:2400)); file = "");
rad = HyperspectralRadiation{FT}(WaveLengthSet{FT}(collect(400:50:2400)); file = ClimaCache.FILE_SUN);
```
"""
HyperspectralRadiation{FT}(wls::WaveLengthSet = WaveLengthSet{FT}(); file::String = FILE_SUN) where {FT<:AbstractFloat} = (
    @unpack SΛ, NΛ = wls;

    # create arrays
    _e_direct  = zeros(FT, NΛ);
    _e_diffuse = zeros(FT, NΛ);

    # Read data from SCOPE MAT file if file is given, if not use 0
    if isfile(file)
        _suni    = matread(file)["sun"];
        _wl      = _suni["wl"      ];
        __e_dir  = _suni["Edirect" ];
        __e_diff = _suni["Ediffuse"];

        # fill in the arrays
        for _i in 1:NΛ
            _wi = findall( SΛ[_i] .<= _wl .< SΛ[_i+1] );
            if length(_wi)==0 @warn "Some wavelengths out of bounds $(string(SΛ[_i]))" end;
            _e_direct[_i]  = mean( __e_dir[_wi]);
            _e_diffuse[_i] = mean(__e_diff[_wi]);
        end;
    end;

    return HyperspectralRadiation{FT}(_e_diffuse, _e_direct)
);
