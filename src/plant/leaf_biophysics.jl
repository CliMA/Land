#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2021-Aug-04: refactor the structure with constants, variables, and temporary cache
#     2021-Aug-04: add concentrations and characteristic curves altogether
#     2021-Aug-10: add CBC and PRO supoort
#     2021-Agu-10: add constructors within the structure rather than initialize it externally
#     2021-Sep-30: rename LeafBio to LeafBiophysics to be more specific
#     2021-Oct-19: sort variable to prognostic and dignostic catergories
#     2021-Oct-21: rename f_sense and K_SENES to brown and K_BROWN
#     2021-Nov-24: tease apart the characteristic absorption curves to HyperspectralAbsorption
#
#######################################################################################################################################################################################################
"""
$(TYPEDEF)

Struct that contains leaf biophysical traits used to run leaf reflection and transmittance.

# Fields
$(TYPEDFIELDS)
"""
mutable struct LeafBiophysics{FT<:AbstractFloat}
    # parameters that do not change with time
    "Leaf mesophyll structural parameter that describes mesophyll reflection and transmittance"
    MESOPHYLL_N::FT
    "Doubling adding layers"
    NDUB::Int

    # prognostic variables that change with time
    "Anthocynanin content `[ug cm⁻²]`"
    ant::FT
    "Senescent material (brown pigments) fraction `[-]`"
    brown::FT
    "Chlorophyll a and b content `[ug cm⁻²]`"
    cab::FT
    "Carotenoid content `[ug cm⁻²]`"
    car::FT
    "Carbon-based constituents in lma `[g cm⁻²]`"
    cbc::FT
    "Zeaxanthin fraction in Carotenoid (1=all Zeaxanthin, 0=all Violaxanthin) `[-]`"
    f_zeax
    "Equivalent water thickness `[cm]`"
    l_H₂O::FT
    "Dry matter content (dry leaf mass per unit area) `[g cm⁻²]`"
    lma::FT
    "Protein content in lma (pro = lma - cbc) `[g cm⁻²]`"
    pro::FT

    # dignostic variables that change with time
    "Specific absorption coefficients of all materials"
    k_all::Vector{FT}
    "Fluorescence excitation matrix backwards `[-]`"
    mat_b::Matrix{FT}
    "Fluorescence excitation matrix forwards `[-]`"
    mat_f::Matrix{FT}
    "Relative absorption by Chlorophyll `[-]`"
    α_cab::Vector{FT}
    "Relative absorption by Chlorophyll+Carotenoid `[-]`"
    α_cabcar::Vector{FT}
    "Shortwave absorption, 1 .- ρ_SW .- τ_SW  `[-]`"
    α_SW::Vector{FT}
    "Broadband thermal reflectance, related to blackbody emittance `[-]`"
    ρ_LW::FT
    "Shortwave leaf reflectance `[-]`"
    ρ_SW::Vector{FT}
    "Broadband thermal transmission, related to blackbody emittance `[-]`"
    τ_LW::FT
    "Shortwave leaf transmission `[-]`"
    τ_SW::Vector{FT}
end


"""
    LeafBiophysics{FT}(wls::WaveLengthSet) where {FT<:AbstractFloat}

Constructor for [`LeafBiophysics`](@ref), given
- `wls` [`WaveLengthSet`](@ref) type structure

---
# Examples
```julia
lbio = LeafBiophysics{FT}();
lbio = LeafBiophysics{FT}(WaveLengthSet{FT}(collect(400:50:2400)));
```
"""
LeafBiophysics{FT}(wls::WaveLengthSet = WaveLengthSet{FT}()) where {FT<:AbstractFloat} = (
    @unpack NΛ, NΛ_SIF, NΛ_SIFE, SΛ = wls;

    return LeafBiophysics{FT}(1.4, 10, 0, 0, 40, 10, 0, 0, 0.01, 0.012, 0, zeros(FT,NΛ), zeros(FT,NΛ_SIF,NΛ_SIFE), zeros(FT,NΛ_SIF,NΛ_SIFE), zeros(FT,NΛ), zeros(FT,NΛ), zeros(FT,NΛ), 0.01,
                              zeros(FT,NΛ), 0.01, zeros(FT,NΛ))
);
