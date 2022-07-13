#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jun-14: add abstract type for soil albedo
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractSoilAlbedo:
- [`BroadbandSoilAlbedo`](@ref)
- [`HyperspectralSoilAlbedo`](@ref)

"""
abstract type AbstractSoilAlbedo{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-14: add struct for broadband soil albedo
#     2022-Jun-14: make soil albedo a two-element vector for PAR and NIR
#     2022-Jul-13: use @kwdef for the constructor
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for broadband soil albedo

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct BroadbandSoilAlbedo{FT<:AbstractFloat} <: AbstractSoilAlbedo{FT}
    # diagnostic variables that change with time
    "Net diffuse radiation at top soil `[W m⁻²]`"
    e_net_diffuse::FT = 0
    "Net direct radiation at top soil `[W m⁻²]`"
    e_net_direct::FT = 0
    "Net longwave energy absorption `[W m⁻²]`"
    r_net_lw::FT = 0
    "Net shortwave energy absorption `[W m⁻²]`"
    r_net_sw::FT = 0
    "Reflectance for longwave radiation"
    ρ_lw::FT = 0.06
    "Reflectance for shortwave radiation (for PAR and NIR)"
    ρ_sw::Vector{FT} = FT[0,0]
end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-14: add struct for hyperspectral soil albedo
#     2022-Jun-14: add constructor
#     2022-Jun-14: add fields to compute soil hyperspectral albedo in CanopyRadiativeTransfer.jl
#     2022-Jun-14: add wls in constructor function and remove n_λ
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for hyperspectral soil albedo

# Fields

$(TYPEDFIELDS)

"""
mutable struct HyperspectralSoilAlbedo{FT<:AbstractFloat} <: AbstractSoilAlbedo{FT}
    # diagnostic variables that change with time
    "Net diffuse radiation at top soil `[mW m⁻² nm⁻¹]`"
    e_net_diffuse::Vector{FT}
    "Net direct radiation at top soil `[mW m⁻² nm⁻¹]`"
    e_net_direct::Vector{FT}
    "A matrix of characteristic curves"
    mat_ρ::Matrix{FT}
    "Net longwave energy absorption `[W m⁻²]`"
    r_net_lw::FT
    "Net shortwave energy absorption `[W m⁻²]`"
    r_net_sw::FT
    "Reflectance for longwave radiation"
    ρ_lw::FT
    "Reflectance for shortwave radiation"
    ρ_sw::Vector{FT}

    # caches to speed up calculations
    "Cache variable with length of NIR"
    _tmp_vec_nir::Vector{FT}
    "Weights of the four characteristic curves"
    _weight::Vector{FT}
    "Cache variable to store ρ_PAR and ρ_NIR (a segmented curve)"
    _ρ_sw::Vector{FT}
end


"""

    HyperspectralSoilAlbedo{FT}(wls::WaveLengthSet{FT}= WaveLengthSet{FT}()) where {FT<:AbstractFloat}

Construct a hyperspectral soil albedo struct, given
- `wls` [`WaveLengthSet`](@ref) type struct that defines wavelength settings

"""
HyperspectralSoilAlbedo{FT}(wls::WaveLengthSet{FT}= WaveLengthSet{FT}()) where {FT<:AbstractFloat} = (
    @unpack IΛ_NIR, NΛ, SΛ = wls;

    # read data that has a 10 nm stepping
    _dat = read_csv(SOIL_GSV);
    _raw_4 = Matrix{FT}(_dat[:,2:end-1]);

    # extend the data to 1 nm stepping by interpolating the data
    _wlr = _dat.WL;
    _wle = collect(400:2500.1);
    _ext_4 = Matrix{FT}(undef, length(_wle), 4);
    for _ie in 1:size(_ext_4,1)-1
        _wl = _wle[_ie];
        _ir = Int(fld(_wl - 400, 10)) + 1;
        _a = (_wl - _wlr[_ir]) * 0.1;
        _ext_4[_ie,:] .= (1-_a) .* _raw_4[_ir,:] .+ _a .* _raw_4[_ir+1,:];
    end;
    _ext_4[end,:] = _raw_4[end,:];

    # rescale the data to match the steppings of wavelength set
    _res_4 = Matrix{FT}(undef, NΛ, 4);

    # fill in the arrays
    for _i_res in 1:NΛ
        _wo = findall( (_wle .>= SΛ[_i_res]) .& (_wle .< SΛ[_i_res+1]) )
        if length(_wo) == 0
            @warn "Some wavelengths out of bounds $(string(SΛ[_i_res]))";
        end;
        _res_4[_i_res,1] = mean( _ext_4[_wo,1] );
        _res_4[_i_res,2] = mean( _ext_4[_wo,2] );
        _res_4[_i_res,3] = mean( _ext_4[_wo,3] );
        _res_4[_i_res,4] = mean( _ext_4[_wo,4] );
    end;

    return HyperspectralSoilAlbedo{FT}(
                zeros(FT,NΛ),               # e_net_diffuse
                zeros(FT,NΛ),               # e_net_direct
                _res_4,                     # mat_ρ
                0,                          # r_net_lw
                0,                          # r_net_sw
                0.06,                       # ρ_lw
                zeros(FT,NΛ),               # ρ_sw
                zeros(FT,length(IΛ_NIR)),   # _tmp_vec_nir
                zeros(FT,4),                # _weight
                zeros(FT,NΛ)                # _ρ_sw
    )
);


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jul-13: add SoilLayer structure
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for soil layer

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SoilLayer{FT<:AbstractFloat}
    # parameters that do not change with time
    "Soil moisture retention curve"
    VC::Union{BrooksCorey{FT}, VanGenuchten{FT}} = VanGenuchten{FT}("Loam")
    "Mean depth"
    Z::FT = -0.5
    "Depth boundaries"
    ZS::Vector{FT} = FT[0,-1]

    # prognostic variables that change with time
    "Temperature"
    t::FT = T_25()
    "Soil water content"
    θ::FT = VC.Θ_SAT
end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-08: add Soil structure
#     2022-Jun-08: add constructor
#     2022-Jun-09: add fields: e_net_diffuse, e_net_direct
#     2022-Jun-10: add fields: r_net_lw, r_net_sw, ρ_lw
#     2022-Jun-14: add abstractized soil albedo
#     2022-Jun-13: use Union instead of Abstract... for type definition
#     2022-Jun-14: add field for soil color class
#     2022-Jun-14: separate the constructor for hyperspectral albedo
#     2022-Jun-14: separate the constructor for broadband albedo
#     2022-Jul-13: move VC, Z, t, and θ to SoilLayer
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for Soil

# Fields

$(TYPEDFIELDS)

"""
mutable struct Soil{FT<:AbstractFloat}
    # parameters that do not change with time
    "Albedo related structure"
    ALBEDO::Union{BroadbandSoilAlbedo{FT}, HyperspectralSoilAlbedo{FT}}
    "Color class as in CLM"
    COLOR::Int
    "Soil layers"
    LAYERS::Vector{SoilLayer{FT}}
end


"""

    Soil{FT}(zs::Vector{FT}, wls::WaveLengthSet{FT} = WaveLengthSet{FT}(); soil_type::String = "Loam") where {FT<:AbstractFloat}

Construct a soil struct with hyperspectral albedo, given
- `zs` Soil upper and lower boundaries
- `wls` [`WaveLengthSet`](@ref) type struct that defines wavelength settings
- `soil_type` Soil type name

"""
Soil{FT}(zs::Vector{FT}, wls::WaveLengthSet{FT} = WaveLengthSet{FT}(); soil_type::String = "Loam") where {FT<:AbstractFloat} = (
    _layers = SoilLayer{FT}[];
    for _i in 1:length(zs-1)
        push!(_layers, SoilLayer{FT}(VC = VanGenuchten{FT}(soil_type), Z = mean(zs[_i:_i+1]), ZS = zs[_i:_i+1]));
    end;

    _sab = HyperspectralSoilAlbedo{FT}(wls);

    return Soil{FT}(_sab, 1, _layers)
);


"""

    Soil{FT}(zs::Vector{FT}, broadband::Bool; soil_type::String = "Loam") where {FT<:AbstractFloat}

Construct a soil struct with broadband albedo, given
- `zs` Soil upper and lower boundaries
- `broadband` Indicating broadband soil albedo
- `soil_type` Soil type name

"""
Soil{FT}(zs::Vector{FT}, broadband::Bool; soil_type::String = "Loam") where {FT<:AbstractFloat} = (
    _layers = SoilLayer{FT}[];
    for _i in 1:length(zs-1)
        push!(_layers, SoilLayer{FT}(VC = VanGenuchten{FT}(soil_type), Z = mean(zs[_i:_i+1]), ZS = zs[_i:_i+1]));
    end;

    _sab = BroadbandSoilAlbedo{FT}();

    return Soil{FT}(_sab, 1, _layers)
);
