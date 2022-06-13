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
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for Verhoef LIDF algorithm

# Fields

$(TYPEDFIELDS)

"""
mutable struct CanopyOpticalProperty{FT<:AbstractFloat}
    # diagnostic variables that change with time
    "Backward diffuse->diffuse scatter weight"
    ddb::FT
    "Forward diffuse->diffuse scatter weight"
    ddf::FT
    "Backward diffuse->observer scatter weight"
    dob::FT
    "Forward diffuse->observer scatter weight"
    dof::FT
    "Conversion factor fo for angle towards observer at different inclination and azimuth angles"
    fo::Matrix{FT}
    "Conversion factor fs for angles from solar at different inclination and azimuth angles"
    fs::Matrix{FT}
    "Observer direction beam extinction coefficient weight (diffuse)"
    ko::FT
    "Solar direction beam extinction coefficient weight (direct)"
    ks::FT
    "Probability of directly viewing a leaf in solar direction at different layers"
    p_sunlit::Vector{FT}
    "Probability of directly viewing a leaf in observer direction at different layer boundaries"
    po::Vector{FT}
    "Probability of directly viewing a leaf in solar direction at different layer boundaries"
    ps::Vector{FT}
    "Bi-directional probability of directly viewing a leaf at different layer boundaries (solar->canopy->observer)"
    pso::Vector{FT}
    "Directional->diffuse backscatter weight"
    sdb::FT
    "Directional->diffuse forward scatter weight"
    sdf::FT
    "Solar directional->observer weight of specular2directional backscatter coefficient"
    sob::FT
    "Solar directional->observer weight of specular2directional forward coefficient"
    sof::FT
    "Effective emissivity for different layers"
    ϵ::Vector{FT}
    "Effective reflectance for diffuse->diffuse"
    ρ_dd::Matrix{FT}
    "Effective reflectance for longwave radiation"
    ρ_lw::Vector{FT}
    "Effective reflectance for directional->diffuse"
    ρ_sd::Matrix{FT}
    "Backward scattering coefficient for diffuse->diffuse at different layers and wavelength bins"
    σ_ddb::Matrix{FT}
    "Forward scattering coefficient for diffuse->diffuse at different layers and wavelength bins"
    σ_ddf::Matrix{FT}
    "Backward scattering coefficient for diffuse->observer at different layers and wavelength bins"
    σ_dob::Matrix{FT}
    "Forward scattering coefficient for diffuse->observer at different layers and wavelength bins"
    σ_dof::Matrix{FT}
    "Backward scattering coefficient for solar directional->diffuse at different layers and wavelength bins"
    σ_sdb::Matrix{FT}
    "Forward scattering coefficient for solar directional->diffuse at different layers and wavelength bins"
    σ_sdf::Matrix{FT}
    "Bidirectional from solar to observer scattering coefficient at different layers and wavelength bins"
    σ_so::Matrix{FT}
    "Effective tranmittance for diffuse->diffuse"
    τ_dd::Matrix{FT}
    "Effective tranmittance for longwave radiation"
    τ_lw::Vector{FT}
    "Effective tranmittance for solar directional->diffuse"
    τ_sd::Matrix{FT}

    # caches to speed up calculations
    "cos(inclination) * cos(vza) at different inclination angles"
    _Co::Vector{FT}
    "cos(inclination) * cos(sza) at different inclination angles"
    _Cs::Vector{FT}
    "sin(inclination) * sin(vza) at different inclination angles"
    _So::Vector{FT}
    "sin(inclination) * sin(sza) at different inclination angles"
    _Ss::Vector{FT}
    "abs of fo"
    _abs_fo::Matrix{FT}
    "abs of fs"
    _abs_fs::Matrix{FT}
    "abs of fs * fo"
    _abs_fs_fo::Matrix{FT}
    "Weighted sum of cos²(inclination)"
    _bf::FT
    "Cosine of Θ_AZI - raa"
    _cos_θ_azi_raa::Vector{FT}
    "fo * cos Θ_INCL"
    _fo_cos_θ_incl::Matrix{FT}
    "fs * cos Θ_INCL"
    _fs_cos_θ_incl::Matrix{FT}
    "fs * fo"
    _fs_fo::Matrix{FT}
    "Outgoing beam extinction coefficient weights at different inclination angles"
    _ko::Vector{FT}
    "Solar beam extinction coefficient weights at different inclination angles"
    _ks::Vector{FT}
    "Upwelling matrix for SIF excitation"
    _mat⁺::Matrix{FT}
    "Downwelling matrix for SIF excitation"
    _mat⁻::Matrix{FT}
    "Backward scattering coefficients at different inclination angles"
    _sb::Vector{FT}
    "Forward scattering coefficients at different inclination angles"
    _sf::Vector{FT}
    "Temporary cache used for matrix adding up purpose (n_incl * n_azi)"
    _tmp_mat_incl_azi_1::Matrix{FT}
    "Temporary cache used for matrix adding up purpose (n_incl * n_azi)"
    _tmp_mat_incl_azi_2::Matrix{FT}
    "Temporary cache used for vector operations (n_azi)"
    _tmp_vec_azi::Vector{FT}
    "Temporary cache used for vector operations (n_layer)"
    _tmp_vec_layer::Vector{FT}
    "Cache variable to store the SIF information"
    _tmp_vec_sif_1::Vector{FT}
    "Cache variable to store the SIF information"
    _tmp_vec_sif_2::Vector{FT}
    "Cache variable to store the SIF information"
    _tmp_vec_sif_3::Vector{FT}
    "Cache variable to store the SIF information"
    _tmp_vec_sif_4::Vector{FT}
    "Cache variable to store the SIF information"
    _tmp_vec_sif_5::Vector{FT}
    "Cache variable to store the SIF information"
    _tmp_vec_sif_6::Vector{FT}
    "Cache variable to store the SIF excitation information"
    _tmp_vec_sife_1::Vector{FT}
    "Cache variable to store the SIF excitation information"
    _tmp_vec_sife_2::Vector{FT}
    "Cache variable to store the SIF excitation information"
    _tmp_vec_sife_3::Vector{FT}
    "Temporary cache used for vector operations (n_λ)"
    _tmp_vec_λ::Vector{FT}
    "Reflectance for diffuse->diffuse at each canopy layer"
    _ρ_dd::Matrix{FT}
    "Reflectance for longwave radiation at each canopy layer"
    _ρ_lw::Vector{FT}
    "Reflectance for solar directional->diffuse at each canopy layer"
    _ρ_sd::Matrix{FT}
    "Tranmittance for diffuse->diffuse at each canopy layer"
    _τ_dd::Matrix{FT}
    "Tranmittance for longwave radiation at each canopy layer"
    _τ_lw::Vector{FT}
    "Tranmittance for solar directional->diffuse at each canopy layer"
    _τ_sd::Matrix{FT}
    "Tranmittance for solar directional->directional at each canopy layer"
    _τ_ss::FT
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-07: add constructor
#     2022-Jun-07: add more fields: fo, fs, po, ps, pso, _Co, _Cs, _So, _Ss, _abs_fo, _abs_fs, _abs_fs_fo, _cos_θ_azi_raa, _fs_fo, _tmp_mat_incl_azi_1, _tmp_mat_incl_azi_2
#     2022-Jun-08: add more fields: ρ_dd, ρ_sd, σ_ddb, σ_ddf, σ_dob, σ_dof, σ_sdb, σ_sdf, σ_so, τ_dd, τ_sd, _tmp_vec_λ, _ρ_dd, _ρ_sd, _τ_dd, _τ_sd
#     2022-Jun-09: fix documentation
#     2022-Jun-09: rename variables to be more descriptive
#     2022-Jun-10: add more fields: p_sunlit, ϵ, ρ_lw, τ_lw, _mat_down, _mat_down, _tmp_vec_azi, _ρ_lw, _τ_lw
#     2022-Jun-10: add more fields for sif calculations
#     2022-Jun-13: add more fields for sif calculations
#     2022-Jun-13: remove unnecessary cache variables
#
#######################################################################################################################################################################################################
"""

    CanopyOpticalProperty{FT}(; n_azi::Int = 36, n_incl::Int = 9, n_layer::Int = 20, n_λ::Int = 114, n_λe::Int = 45, n_λf::Int = 29) where {FT<:AbstractFloat}

Construct a struct to store canopy optical properties, given
- `n_azi` Number of azimuth angles
- `n_incl` Number of inclination angles
- `n_layer` Number of canopy layers
- `n_λ` Number of wavelength bins
- `n_λe` Number of SIF excitation wavelength bins
- `n_λf` Number of SIF wavelength bins
"""
CanopyOpticalProperty{FT}(; n_azi::Int = 36, n_incl::Int = 9, n_layer::Int = 20, n_λ::Int = 114, n_λe::Int = 45, n_λf::Int = 29) where {FT<:AbstractFloat} = (
    return CanopyOpticalProperty{FT}(
                0,                          # ddb
                0,                          # ddf
                0,                          # dob
                0,                          # dof
                zeros(FT,n_incl,n_azi),     # fo
                zeros(FT,n_incl,n_azi),     # fs
                0,                          # ko
                0,                          # ks
                zeros(FT,n_layer),          # p_sunlit
                zeros(FT,n_layer+1),        # po
                zeros(FT,n_layer+1),        # ps
                zeros(FT,n_layer+1),        # pso
                0,                          # sdb
                0,                          # sdf
                0,                          # sob
                0,                          # sof
                zeros(FT,n_layer),          # ϵ
                zeros(FT,n_λ,n_layer+1),    # ρ_dd
                zeros(FT,n_layer+1),        # ρ_lw
                zeros(FT,n_λ,n_layer+1),    # ρ_sd
                zeros(FT,n_λ,n_layer),      # σ_ddb
                zeros(FT,n_λ,n_layer),      # σ_ddf
                zeros(FT,n_λ,n_layer),      # σ_dob
                zeros(FT,n_λ,n_layer),      # σ_dof
                zeros(FT,n_λ,n_layer),      # σ_sdb
                zeros(FT,n_λ,n_layer),      # σ_sdf
                zeros(FT,n_λ,n_layer),      # σ_so
                zeros(FT,n_λ,n_layer),      # τ_dd
                zeros(FT,n_layer),          # τ_lw
                zeros(FT,n_λ,n_layer),      # τ_sd
                zeros(FT,n_incl),           # _Co
                zeros(FT,n_incl),           # _Cs
                zeros(FT,n_incl),           # _So
                zeros(FT,n_incl),           # _Ss
                zeros(FT,n_incl,n_azi),     # _abs_fo
                zeros(FT,n_incl,n_azi),     # _abs_fs
                zeros(FT,n_incl,n_azi),     # _abs_fs_fo
                0,                          # _bf
                zeros(FT,n_azi),            # _cos_θ_azi_raa
                zeros(FT,n_incl,n_azi),     # _fo_cos_θ_incl
                zeros(FT,n_incl,n_azi),     # _fs_cos_θ_incl
                zeros(FT,n_incl,n_azi),     # _fs_fo
                zeros(FT,n_incl),           # _ko
                zeros(FT,n_incl),           # _ks
                zeros(FT,n_λf,n_λe),        # _mat⁺
                zeros(FT,n_λf,n_λe),        # _mat⁻
                zeros(FT,n_incl),           # _sb
                zeros(FT,n_incl),           # _sf
                zeros(FT,n_incl,n_azi),     # _tmp_mat_incl_azi_1
                zeros(FT,n_incl,n_azi),     # _tmp_mat_incl_azi_2
                zeros(FT,n_azi),            # _tmp_vec_azi
                zeros(FT,n_layer),          # _tmp_vec_layer
                zeros(FT,n_λf),             # _tmp_vec_sif_1
                zeros(FT,n_λf),             # _tmp_vec_sif_2
                zeros(FT,n_λf),             # _tmp_vec_sif_3
                zeros(FT,n_λf),             # _tmp_vec_sif_4
                zeros(FT,n_λf),             # _tmp_vec_sif_5
                zeros(FT,n_λf),             # _tmp_vec_sif_6
                zeros(FT,n_λe),             # _tmp_vec_sife_1
                zeros(FT,n_λe),             # _tmp_vec_sife_2
                zeros(FT,n_λe),             # _tmp_vec_sife_3
                zeros(FT,n_λ),              # _tmp_vec_λ
                zeros(FT,n_λ,n_layer),      # _ρ_dd
                zeros(FT,n_layer),          # _ρ_lw
                zeros(FT,n_λ,n_layer),      # _ρ_sd
                zeros(FT,n_λ,n_layer),      # _τ_dd
                zeros(FT,n_layer),          # _τ_lw
                zeros(FT,n_λ,n_layer),      # _τ_sd
                0                           # _τ_ss
    )
);
