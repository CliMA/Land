###############################################################################
#
# Leaf optical parameters
#
###############################################################################
"""
    struct LeafOpticals{FT}

Struct for leaf optical properties using Array

# Fields
$(DocStringExtensions.FIELDS)
"""
struct LeafOpticals{FT}
    # TODO Add explanations to each field
    nr    ::Array{FT,1}
    Km    ::Array{FT,1}
    Kab   ::Array{FT,1}
    Kant  ::Array{FT,1}
    Kcar  ::Array{FT,1}
    Kw    ::Array{FT,1}
    KBrown::Array{FT,1}
    phi   ::Array{FT,1}
    KcaV  ::Array{FT,1}
    KcaZ  ::Array{FT,1}
    lambda::Array{FT,1}
end




"""
    create_leaf_opticals(swl::Array, file::String)

Create an `AbstractLeafOptiPara` struct, given
- `swl` Standard wave length
- `file` Input file name
"""
function create_leaf_opticals(swl::Array{FT,1}, file::String=file_Opti) where {FT<:AbstractFloat}
    N = length(swl)-1

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
    # println("Reading Optical Parameters from ", swl[1], " to ", swl[end], " length: ", length(swl))
    @inbounds for i in 1:N
        wo = findall( (_lambda.>=swl[i]) .& (_lambda.<swl[i+1]) )
        if length(wo)==0
            println("Warning, some wavelengths out of bounds ", swl[i])
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
    return LeafOpticals{FT}(nr, Km, Kab, Kant, Kcar, Kw, KBrown, phi, KcaV, KcaZ, lambda)
end
