#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-07: add CanopyOptics struct (will be a field for canopy structure)
#     2022-Jun-07: add more fields: fo, fs, po, ps, pso, _Co, _Cs, _So, _Ss, _abs_fo, _abs_fs, _abs_fs_fo, _cos_θ_azi_raa, _fs_fo, _tmp_mat_incl_azi_1, _tmp_mat_incl_azi_2
#     2022-Jun-08: add more fields: ρ_dd, ρ_sd, σ_ddb, σ_ddf, σ_dob, σ_dof, σ_sdb, σ_sdf, σ_so, τ_dd, τ_sd, _tmp_vec_λ, _ρ_dd, _ρ_sd, _τ_dd, _τ_sd
#     2022-Jun-09: rename variables to be more descriptive
#     2022-Jun-10: add more fields: p_sunlit, ϵ, ρ_lw, τ_lw, _mat_down, _mat_down, _tmp_vec_azi, _ρ_lw, _τ_lw
#     2022-Jun-10: add more fields for sif calculations
#     2022-Jun-13: add more fields for sif calculations
#     2022-Jun-13: remove unnecessary cache variables
#     2022-Jun-15: rename to HyperspectralMLCanopyOpticalProperty
#     2022-Jul-19: use kwdef for the constructor
#     2022-Jul-19: add dimension control to struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for Verhoef LIDF algorithm

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct HyperspectralMLCanopyOpticalProperty{FT<:AbstractFloat}
    # dimensions
    "Dimension of azimuth angles"
    DIM_AZI::Int = 36
    "Dimension of inclination angles"
    DIM_INCL::Int = 9
    "Dimension of canopy layers"
    DIM_LAYER::Int = 20
    "Dimension of SIF wave length bins"
    DIM_SIF::Int = 29
    "Dimension of SIF excitation wave length bins"
    DIM_SIFE::Int = 45
    "Dimension of short wave length bins"
    DIM_WL::Int = 114

    # diagnostic variables that change with time
    "Backward diffuse->diffuse scatter weight"
    ddb::FT = 0
    "Forward diffuse->diffuse scatter weight"
    ddf::FT = 0
    "Backward diffuse->observer scatter weight"
    dob::FT = 0
    "Forward diffuse->observer scatter weight"
    dof::FT = 0
    "Conversion factor fo for angle towards observer at different inclination and azimuth angles"
    fo::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
    "Conversion factor fs for angles from solar at different inclination and azimuth angles"
    fs::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
    "Observer direction beam extinction coefficient weight (diffuse)"
    ko::FT = 0
    "Solar direction beam extinction coefficient weight (direct)"
    ks::FT = 0
    "Probability of directly viewing a leaf in solar direction at different layers"
    p_sunlit::Vector{FT} = zeros(FT, DIM_LAYER)
    "Probability of directly viewing a leaf in observer direction at different layer boundaries"
    po::Vector{FT} = zeros(FT, DIM_LAYER+1)
    "Probability of directly viewing a leaf in solar direction at different layer boundaries"
    ps::Vector{FT} = zeros(FT, DIM_LAYER+1)
    "Bi-directional probability of directly viewing a leaf at different layer boundaries (solar->canopy->observer)"
    pso::Vector{FT} = zeros(FT, DIM_LAYER+1)
    "Directional->diffuse backscatter weight"
    sdb::FT = 0
    "Directional->diffuse forward scatter weight"
    sdf::FT = 0
    "Solar directional->observer weight of specular2directional backscatter coefficient"
    sob::FT = 0
    "Solar directional->observer weight of specular2directional forward coefficient"
    sof::FT = 0
    "Effective emissivity for different layers"
    ϵ::Vector{FT} = zeros(FT, DIM_LAYER)
    "Effective reflectance for diffuse->diffuse"
    ρ_dd::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER+1)
    "Effective reflectance for longwave radiation"
    ρ_lw::Vector{FT} = zeros(FT, DIM_LAYER+1)
    "Effective reflectance for directional->diffuse"
    ρ_sd::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER+1)
    "Backward scattering coefficient for diffuse->diffuse at different layers and wavelength bins"
    σ_ddb::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Forward scattering coefficient for diffuse->diffuse at different layers and wavelength bins"
    σ_ddf::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Backward scattering coefficient for diffuse->observer at different layers and wavelength bins"
    σ_dob::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Forward scattering coefficient for diffuse->observer at different layers and wavelength bins"
    σ_dof::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Backward scattering coefficient for solar directional->diffuse at different layers and wavelength bins"
    σ_sdb::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Forward scattering coefficient for solar directional->diffuse at different layers and wavelength bins"
    σ_sdf::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Bidirectional from solar to observer scattering coefficient at different layers and wavelength bins"
    σ_so::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Effective tranmittance for diffuse->diffuse"
    τ_dd::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Effective tranmittance for longwave radiation"
    τ_lw::Vector{FT} = zeros(FT, DIM_LAYER)
    "Effective tranmittance for solar directional->diffuse"
    τ_sd::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)

    # caches to speed up calculations
    "cos(inclination) * cos(vza) at different inclination angles"
    _Co::Vector{FT} = zeros(FT, DIM_INCL)
    "cos(inclination) * cos(sza) at different inclination angles"
    _Cs::Vector{FT} = zeros(FT, DIM_INCL)
    "sin(inclination) * sin(vza) at different inclination angles"
    _So::Vector{FT} = zeros(FT, DIM_INCL)
    "sin(inclination) * sin(sza) at different inclination angles"
    _Ss::Vector{FT} = zeros(FT, DIM_INCL)
    "abs of fo"
    _abs_fo::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
    "abs of fs"
    _abs_fs::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
    "abs of fs * fo"
    _abs_fs_fo::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
    "Weighted sum of cos²(inclination)"
    _bf::FT = 0
    "Cosine of Θ_AZI - raa"
    _cos_θ_azi_raa::Vector{FT} = zeros(FT, DIM_AZI)
    "fo * cos Θ_INCL"
    _fo_cos_θ_incl::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
    "fs * cos Θ_INCL"
    _fs_cos_θ_incl::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
    "fs * fo"
    _fs_fo::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
    "Outgoing beam extinction coefficient weights at different inclination angles"
    _ko::Vector{FT} = zeros(FT, DIM_INCL)
    "Solar beam extinction coefficient weights at different inclination angles"
    _ks::Vector{FT} = zeros(FT, DIM_INCL)
    "Upwelling matrix for SIF excitation"
    _mat⁺::Matrix{FT} = zeros(FT, DIM_SIF, DIM_SIFE)
    "Downwelling matrix for SIF excitation"
    _mat⁻::Matrix{FT} = zeros(FT, DIM_SIF, DIM_SIFE)
    "Backward scattering coefficients at different inclination angles"
    _sb::Vector{FT} = zeros(FT, DIM_INCL)
    "Forward scattering coefficients at different inclination angles"
    _sf::Vector{FT} = zeros(FT, DIM_INCL)
    "Temporary cache used for matrix adding up purpose (DIM_INCL * DIM_AZI)"
    _tmp_mat_incl_azi_1::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
    "Temporary cache used for matrix adding up purpose (DIM_INCL * DIM_AZI)"
    _tmp_mat_incl_azi_2::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
    "Temporary cache used for vector operations (DIM_AZI)"
    _tmp_vec_azi::Vector{FT} = zeros(FT, DIM_AZI)
    "Temporary cache used for vector operations (DIM_LAYER)"
    _tmp_vec_layer::Vector{FT} = zeros(FT, DIM_LAYER)
    "Cache variable to store the SIF information"
    _tmp_vec_sif_1::Vector{FT} = zeros(FT, DIM_SIF)
    "Cache variable to store the SIF information"
    _tmp_vec_sif_2::Vector{FT} = zeros(FT, DIM_SIF)
    "Cache variable to store the SIF information"
    _tmp_vec_sif_3::Vector{FT} = zeros(FT, DIM_SIF)
    "Cache variable to store the SIF information"
    _tmp_vec_sif_4::Vector{FT} = zeros(FT, DIM_SIF)
    "Cache variable to store the SIF information"
    _tmp_vec_sif_5::Vector{FT} = zeros(FT, DIM_SIF)
    "Cache variable to store the SIF information"
    _tmp_vec_sif_6::Vector{FT} = zeros(FT, DIM_SIF)
    "Cache variable to store the SIF excitation information"
    _tmp_vec_sife_1::Vector{FT} = zeros(FT, DIM_SIFE)
    "Cache variable to store the SIF excitation information"
    _tmp_vec_sife_2::Vector{FT} = zeros(FT, DIM_SIFE)
    "Cache variable to store the SIF excitation information"
    _tmp_vec_sife_3::Vector{FT} = zeros(FT, DIM_SIFE)
    "Temporary cache used for vector operations (DIM_WL)"
    _tmp_vec_λ::Vector{FT} = zeros(FT, DIM_WL)
    "Reflectance for diffuse->diffuse at each canopy layer"
    _ρ_dd::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Reflectance for longwave radiation at each canopy layer"
    _ρ_lw::Vector{FT} = zeros(FT, DIM_LAYER)
    "Reflectance for solar directional->diffuse at each canopy layer"
    _ρ_sd::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Tranmittance for diffuse->diffuse at each canopy layer"
    _τ_dd::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Tranmittance for longwave radiation at each canopy layer"
    _τ_lw::Vector{FT} = zeros(FT, DIM_LAYER)
    "Tranmittance for solar directional->diffuse at each canopy layer"
    _τ_sd::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Tranmittance for solar directional->directional at each canopy layer"
    _τ_ss::FT = 0
end
