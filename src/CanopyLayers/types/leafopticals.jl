###############################################################################
#
# Leaf optical parameters
#
###############################################################################
"""
    mutable struct LeafOpticals{FT}

Struct for leaf optical properties

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct LeafOpticals{FT}
    # TODO Add explanations to each field
    nr::Vector{FT} = zeros(FT, 2)
    Km::Vector{FT} = zeros(FT, 2)
    Kab::Vector{FT} = zeros(FT, 2)
    Kant::Vector{FT} = zeros(FT, 2)
    Kcar::Vector{FT} = zeros(FT, 2)
    Kw::Vector{FT} = zeros(FT, 2)
    KBrown::Vector{FT} = zeros(FT, 2)
    phi::Vector{FT} = zeros(FT, 2)
    KcaV::Vector{FT} = zeros(FT, 2)
    KcaZ::Vector{FT} = zeros(FT, 2)
    "Wave length `[nm]`, same as `WL` in [`WaveLengths`](@ref)`"
    lambda::Vector{FT} = zeros(FT, 2)
end

LeafOpticals{FT}(sWL::Vector{FT}, opfn::String) where {FT} = (
    N = length(sWL) - 1;

    # reading data
    _opti   = matread(opfn)["optipar"];
    _nr     = _opti["nr"  ];
    _Km     = _opti["Kdm" ];
    _Kab    = _opti["Kab" ];
    _Kant   = _opti["Kant"];
    _Kcar   = _opti["Kca" ];
    _Kw     = _opti["Kw"  ];
    _KBrown = _opti["Ks"  ];
    _phi    = _opti["phi" ];
    _KcaV   = _opti["KcaV"];
    _KcaZ   = _opti["KcaZ"];
    _lambda = _opti["wl"  ];

    # create data to parse
    nr     = zeros(FT, N);
    Km     = zeros(FT, N);
    Kab    = zeros(FT, N);
    Kant   = zeros(FT, N);
    Kcar   = zeros(FT, N);
    Kw     = zeros(FT, N);
    KBrown = zeros(FT, N);
    phi    = zeros(FT, N);
    KcaV   = zeros(FT, N);
    KcaZ   = zeros(FT, N);
    lambda = zeros(FT, N);

    # fill in the data arrays
    # println("Reading Optical Parameters from ", sWL[1], " to ", sWL[end])
    @inbounds for i in 1:N
        wo = findall( (_lambda.>=sWL[i]) .& (_lambda.<sWL[i+1]) );
        if length(wo)==0
            @warn "Warning, some wavelengths out of bounds $(string(sWL[i]))";
        end

        nr[i]     = mean(    _nr[wo]);
        Km[i]     = mean(    _Km[wo]);
        Kab[i]    = mean(   _Kab[wo]);
        Kant[i]   = mean(  _Kant[wo]);
        Kcar[i]   = mean(  _Kcar[wo]);
        Kw[i]     = mean(    _Kw[wo]);
        KBrown[i] = mean(_KBrown[wo]);
        phi[i]    = mean(   _phi[wo]);
        KcaV[i]   = mean(  _KcaV[wo]);
        KcaZ[i]   = mean(  _KcaZ[wo]);
        lambda[i] = mean(_lambda[wo]);
    end;

    # return the created struct
    return LeafOpticals{FT}(nr, Km, Kab, Kant, Kcar, Kw, KBrown, phi, KcaV, KcaZ, lambda)
);
