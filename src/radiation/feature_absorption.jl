#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2021-Aug-04: refactor the structure with constants, variables, and temporary cache
#     2021-Aug-04: add concentrations and characteristic curves altogether
#     2021-Aug-10: add CBC and PRO supoort
#     2021-Agu-10: add constructors within the structure rather than initialize it externally
#     2021-Sep-30: rename LeafBio to LeafBiophysics to be more specific
#     2021-Oct-19: sort variable to prognostic and dignostic catergories
#     2021-Oct-21: rename f_sense and K_SENES to brown and K_BROWN
#     2021-Nov-24: tease apart the characteristic absorption curves to HyperspectralAbsorption
#     2022-May-25: fix documentation
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Immutable struct that contains leaf biophysical traits used to run leaf reflection and transmittance.

# Fields

$(TYPEDFIELDS)

"""
struct HyperspectralAbsorption{FT<:AbstractFloat}
    # parameters that do not change with time
    "Specific absorption coefficients of anthocynanin `[-]`"
    K_ANT::Vector{FT}
    "Specific absorption coefficients of senescent material (brown pigments) `[-]`"
    K_BROWN::Vector{FT}
    "Specific absorption coefficients of chlorophyll a and b `[-]`"
    K_CAB::Vector{FT}
    "Specific absorption coefficients of violaxanthin carotenoid `[-]`"
    K_CAR_V::Vector{FT}
    "Specific absorption coefficients of zeaxanthin carotenoid `[-]`"
    K_CAR_Z::Vector{FT}
    "Specific absorption coefficients of carbon-based constituents `[-]`"
    K_CBC::Vector{FT}
    "Specific absorption coefficients of water `[-]`"
    K_H₂O::Vector{FT}
    "Specific absorption coefficients of dry matter `[-]`"
    K_LMA::Vector{FT}
    "Specific absorption coefficients of protein `[-]`"
    K_PRO::Vector{FT}
    "Specific absorption coefficients of PS I and II `[-]`"
    K_PS::Vector{FT}
    "Refractive index `[-]`"
    NR::Vector{FT}
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2021-Nov-24: add constructor
# To do
#     TODO: use Netcdf for the MAT artifact (as MAT is not easily accessible)
#
#######################################################################################################################################################################################################
"""

    HyperspectralAbsorption{FT}(wls::WaveLengthSet = WaveLengthSet{FT}(); opti::String = OPTI_2021) where {FT<:AbstractFloat}

Constructor for [`HyperspectralAbsorption`](@ref), given
- `wls` [`WaveLengthSet`](@ref) type structure
- `opti` Path to leaf optical properties

"""
HyperspectralAbsorption{FT}(wls::WaveLengthSet = WaveLengthSet{FT}(); opti::String = OPTI_2021) where {FT<:AbstractFloat} = (
    @unpack NΛ, NΛ_SIF, NΛ_SIFE, SΛ = wls;

    # read data from the MAT file
    _opti    = matread(opti)["optipar"];
    __nr     = _opti["nr"  ];
    __Klma   = _opti["Kdm" ];
    __Kcab   = _opti["Kab" ];
    __Kant   = _opti["Kant"];
    __Kh2o   = _opti["Kw"  ];
    __Kbrown = _opti["Ks"  ];
    __Kps    = _opti["phi" ];
    __KcaV   = _opti["KcaV"];
    __KcaZ   = _opti["KcaZ"];
    __lambda = _opti["wl"  ];

    # read protein and carbon-based constituents if exist (not in file OPTI_2017)
    if opti == OPTI_2017
        __Kpro = __Klma;
        __Kcbc = __Klma;
    else
        __Kpro = _opti["Kp"  ];
        __Kcbc = _opti["Kcbc"];
    end;

    # create data to parse
    _nr     = zeros(FT, NΛ);
    _Klma   = zeros(FT, NΛ);
    _Kcab   = zeros(FT, NΛ);
    _Kant   = zeros(FT, NΛ);
    _Kh2o   = zeros(FT, NΛ);
    _Kbrown = zeros(FT, NΛ);
    _Kps    = zeros(FT, NΛ);
    _KcaV   = zeros(FT, NΛ);
    _KcaZ   = zeros(FT, NΛ);
    _Kpro   = zeros(FT, NΛ);
    _Kcbc   = zeros(FT, NΛ);

    # fill in the data arrays
    @inbounds for _i in 1:NΛ
        _wo = findall( SΛ[_i] .<= __lambda .< SΛ[_i+1] );
        if length(_wo) < 1 @warn "Warning, some wavelengths out of bounds $(string(SΛ[_i]))" end;

        _nr[_i]     = mean(__nr[_wo]    );
        _Klma[_i]   = mean(__Klma[_wo]  );
        _Kcab[_i]   = mean(__Kcab[_wo]  );
        _Kant[_i]   = mean(__Kant[_wo]  );
        _Kh2o[_i]   = mean(__Kh2o[_wo]  );
        _Kbrown[_i] = mean(__Kbrown[_wo]);
        _Kps[_i]    = mean(__Kps[_wo]   );
        _KcaV[_i]   = mean(__KcaV[_wo]  );
        _KcaZ[_i]   = mean(__KcaZ[_wo]  );
        _Kpro[_i]   = mean(__Kpro[_wo]  );
        _Kcbc[_i]   = mean(__Kcbc[_wo]  );
    end;

    return HyperspectralAbsorption{FT}(
                _Kant,      # K_ANT
                _Kbrown,    # K_BROWN
                _Kcab,      # K_CAB
                _KcaV,      # K_CAR_V
                _KcaZ,      # K_CAR_Z
                _Kcbc,      # K_CBC
                _Kh2o,      # K_H₂O
                _Klma,      # K_LMA
                _Kpro,      # K_PRO
                _Kps,       # K_PS
                _nr         # NR
    )
);
