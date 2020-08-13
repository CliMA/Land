###############################################################################
#
# Create LeafOpticals
#
###############################################################################
"""
    create_leaf_opticals(sWL::Array, file::String)

Create an `AbstractLeafOptiPara` struct, given
- `sWL` Standard wave length
- `file` Input file name
"""
function create_leaf_opticals(sWL::Array{FT,1}, file::String=file_Opti) where {FT<:AbstractFloat}
    N = length(sWL)-1

    # reading data
    _opti   = matread(file)["optipar"]
    _nr     = _opti["nr"  ]
    _Km     = _opti["Kdm" ]
    _Kab    = _opti["Kab" ]
    _Kant   = _opti["Kant"]
    _Kcar   = _opti["Kca" ]
    _Kw     = _opti["Kw"  ]
    _KBrown = _opti["Ks"  ]
    _phi    = _opti["phi" ]
    _KcaV   = _opti["KcaV"]
    _KcaZ   = _opti["KcaZ"]
    _lambda = _opti["wl"  ]

    # create data to parse
    nr     = zeros(FT, N)
    Km     = zeros(FT, N)
    Kab    = zeros(FT, N)
    Kant   = zeros(FT, N)
    Kcar   = zeros(FT, N)
    Kw     = zeros(FT, N)
    KBrown = zeros(FT, N)
    phi    = zeros(FT, N)
    KcaV   = zeros(FT, N)
    KcaZ   = zeros(FT, N)
    lambda = zeros(FT, N)

    # fill in the data arrays
    # println("Reading Optical Parameters from ", sWL[1], " to ", sWL[end], " length: ", length(sWL))
    @inbounds for i in 1:N
        wo = findall( (_lambda.>=sWL[i]) .& (_lambda.<sWL[i+1]) )
        if length(wo)==0
            println("Warning, some wavelengths out of bounds ", sWL[i])
        end

        nr[i]     = mean(    _nr[wo])
        Km[i]     = mean(    _Km[wo])
        Kab[i]    = mean(   _Kab[wo])
        Kant[i]   = mean(  _Kant[wo])
        Kcar[i]   = mean(  _Kcar[wo])
        Kw[i]     = mean(    _Kw[wo])
        KBrown[i] = mean(_KBrown[wo])
        phi[i]    = mean(   _phi[wo])
        KcaV[i]   = mean(  _KcaV[wo])
        KcaZ[i]   = mean(  _KcaZ[wo])
        lambda[i] = mean(_lambda[wo])
    end

    # return the created struct
    return LeafOpticals{FT}(nr,
                            Km,
                            Kab,
                            Kant,
                            Kcar,
                            Kw,
                            KBrown,
                            phi,
                            KcaV,
                            KcaZ,
                            lambda)
end
