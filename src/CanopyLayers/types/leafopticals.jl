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
    opti_file::String

    # TODO Add explanations to each field
    nr::Vector{FT} = read_nc(opti_file, "NR")
    Km::Vector{FT} = read_nc(opti_file, "K_LMA")
    Kab::Vector{FT} = read_nc(opti_file, "K_CAB")
    Kant::Vector{FT} = read_nc(opti_file, "K_ANT")
    Kw::Vector{FT} = read_nc(opti_file, "K_H₂O")
    KBrown::Vector{FT} = read_nc(opti_file, "K_BROWN")
    phi::Vector{FT} = read_nc(opti_file, "K_PS")
    KcaV::Vector{FT} = read_nc(opti_file, "K_CAR_V")
    KcaZ::Vector{FT} = read_nc(opti_file, "K_CAR_Z")
    Kcar::Vector{FT} = (KcaV .+ KcaZ) ./ 2
end
