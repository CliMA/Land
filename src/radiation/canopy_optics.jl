#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-07: add CanopyOptics struct (will be a field for canopy structure)
#     2022-Jun-07: add more cache fields: fo, fs, po, ps, pso, _Co, _Cs, _So, _Ss, _abs_fo, _abs_fs, _abs_fs_fo, _cos_θ_azi_raa, _fs_fo, _tmp_mat_incl_azi_1, _tmp_mat_incl_azi_2
#     2022-Jun-08: add more cache fields: ρ_dd, ρ_dv, σ_ddb, σ_ddf, σ_dvb, σ_dvf, σ_vdb, σ_vdf, σ_vv, τ_dd, τ_dv, _tmp_vec_λ, _ρ_dd, _ρ_dv, _τ_dd, _τ_dv
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
    "Diffuse -> Diffuse backscatter weight"
    ddb::FT
    "Diffuse -> Diffuse forward scatter weight"
    ddf::FT
    "Diffuse -> Outgoing backscatter weight"
    dob::FT
    "Diffuse -> Outgoing forward scatter weight"
    dof::FT
    "Conversion factor fo for angle towards observer at different inclination and azimuth angles (not sun like fs)"
    fo::Matrix{FT}
    "Conversion factor fs to compute irradiance on inclined leaf at different inclination and azimuth angles"
    fs::Matrix{FT}
    "Outgoing beam extinction coefficient weight"
    ko ::FT
    "Solar beam extinction coefficient weight"
    ks ::FT
    "Probability of directly viewing a leaf in solar direction at different layer boundaries"
    po::Vector{FT}
    "Probability of directly viewing a leaf in viewing direction at different layer boundaries"
    ps::Vector{FT}
    "Bi-directional probability of directly viewing a leaf at different layer boundaries (solar->canopy->viewing)"
    pso::Vector{FT}
    "Solar -> Diffuse backscatter weight"
    sdb::FT
    "Solar -> Diffuse forward scatter weight"
    sdf::FT
    "Solar -> Outgoing weight of specular2directional backscatter coefficient"
    sob::FT
    "Solar -> Outgoing weight of specular2directional forward coefficient"
    sof::FT
    "Effective reflectance for diffuse->diffuse"
    ρ_dd::Matrix{FT}
    "Effective reflectance for diffuse->directional"
    ρ_dv::Matrix{FT}
    "Backward scattering coefficient for diffuse->diffuse at different layers and wavelength bins"
    σ_ddb::Matrix{FT}
    "Forward scattering coefficient for diffuse->diffuse at different layers and wavelength bins"
    σ_ddf::Matrix{FT}
    "Backward scattering coefficient for diffuse->directional at different layers and wavelength bins"
    σ_dvb::Matrix{FT}
    "Forward scattering coefficient for diffuse->directional at different layers and wavelength bins"
    σ_dvf::Matrix{FT}
    "Backward scattering coefficient for directional->diffuse at different layers and wavelength bins"
    σ_vdb::Matrix{FT}
    "Forward scattering coefficient for directional->diffuse at different layers and wavelength bins"
    σ_vdf::Matrix{FT}
    "Bidirectional scattering coefficient at different layers and wavelength bins"
    σ_vv::Matrix{FT}
    "Effective tranmittance for diffuse->diffuse"
    τ_dd::Matrix{FT}
    "Effective tranmittance for diffuse->directional"
    τ_dv::Matrix{FT}

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
    "fs * fo"
    _fs_fo::Matrix{FT}
    "Outgoing beam extinction coefficient weights at different inclination angles"
    _ko::Vector{FT}
    "Solar beam extinction coefficient weights at different inclination angles"
    _ks::Vector{FT}
    "Backward scattering coefficients at different inclination angles"
    _sb::Vector{FT}
    "Forward scattering coefficients at different inclination angles"
    _sf::Vector{FT}
    "Temporary cache used for matrix adding up purpose (n_incl * n_azi)"
    _tmp_mat_incl_azi_1::Matrix{FT}
    "Temporary cache used for matrix adding up purpose (n_incl * n_azi)"
    _tmp_mat_incl_azi_2::Matrix{FT}
    "Temporary cache used for vector operations (n_λ)"
    _tmp_vec_λ::Vector{FT}
    "Reflectance for diffuse->diffuse at each canopy layer"
    _ρ_dd::Matrix{FT}
    "Reflectance for diffuse->directional at each canopy layer"
    _ρ_dv::Matrix{FT}
    "Tranmittance for diffuse->diffuse at each canopy layer"
    _τ_dd::Matrix{FT}
    "Tranmittance for diffuse->directional at each canopy layer"
    _τ_dv::Matrix{FT}
    "Tranmittance for directional->directional at each canopy layer"
    _τ_vv::FT
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-07: add constructor
#     2022-Jun-07: add more cache fields: fo, fs, po, ps, pso, _Co, _Cs, _So, _Ss, _abs_fo, _abs_fs, _abs_fs_fo, _cos_θ_azi_raa, _fs_fo, _tmp_mat_incl_azi_1, _tmp_mat_incl_azi_2
#     2022-Jun-08: add more cache fields: ρ_dd, ρ_dv, σ_ddb, σ_ddf, σ_dvb, σ_dvf, σ_vdb, σ_vdf, σ_vv, τ_dd, τ_dv, _tmp_vec_λ, _ρ_dd, _ρ_dv, _τ_dd, _τ_dv
#
#######################################################################################################################################################################################################
"""

    CanopyOpticalProperty{FT}(; n_azi::Int = 36, n_incl::Int = 9, n_layer::Int = 20, n_λ::Int = 114) where {FT<:AbstractFloat}

Construct a struct to store canopy optical properties
- `n_azi` Number of azimuth angles
- `n_incl` Number of inclination angles
- `n_layer` Number of canopy layers
- `n_λ` Number of wavelength bins
"""
CanopyOpticalProperty{FT}(; n_azi::Int = 36, n_incl::Int = 9, n_layer::Int = 20, n_λ::Int = 114) where {FT<:AbstractFloat} = (
    return CanopyOpticalProperty{FT}(
                0,                          # ddb
                0,                          # ddf
                0,                          # dob
                0,                          # dof
                zeros(FT,n_incl,n_azi),     # fo
                zeros(FT,n_incl,n_azi),     # fs
                0,                          # ko
                0,                          # ks
                zeros(FT,n_layer+1),        # po
                zeros(FT,n_layer+1),        # ps
                zeros(FT,n_layer+1),        # pso
                0,                          # sdb
                0,                          # sdf
                0,                          # sob
                0,                          # sof
                zeros(FT,n_λ,n_layer+1),    # ρ_dd
                zeros(FT,n_λ,n_layer+1),    # ρ_dv
                zeros(FT,n_λ,n_layer),      # σ_ddb
                zeros(FT,n_λ,n_layer),      # σ_ddf
                zeros(FT,n_λ,n_layer),      # σ_dvb
                zeros(FT,n_λ,n_layer),      # σ_dvf
                zeros(FT,n_λ,n_layer),      # σ_vdb
                zeros(FT,n_λ,n_layer),      # σ_vdf
                zeros(FT,n_λ,n_layer),      # σ_vv
                zeros(FT,n_λ,n_layer),      # τ_dd
                zeros(FT,n_λ,n_layer),      # τ_dv
                zeros(FT,n_incl),           # _Co
                zeros(FT,n_incl),           # _Cs
                zeros(FT,n_incl),           # _So
                zeros(FT,n_incl),           # _Ss
                zeros(FT,n_incl,n_azi),     # _abs_fo
                zeros(FT,n_incl,n_azi),     # _abs_fs
                zeros(FT,n_incl,n_azi),     # _abs_fs_fo
                0,                          # _bf
                zeros(FT,n_azi),            # _cos_θ_azi_raa
                zeros(FT,n_incl,n_azi),     # _fs_fo
                zeros(FT,n_incl),           # _ko
                zeros(FT,n_incl),           # _ks
                zeros(FT,n_incl),           # _sb
                zeros(FT,n_incl),           # _sf
                zeros(FT,n_incl,n_azi),     # _tmp_mat_incl_azi_1
                zeros(FT,n_incl,n_azi),     # _tmp_mat_incl_azi_2
                zeros(FT,n_λ),              # _tmp_vec_λ
                zeros(FT,n_λ,n_layer),      # _ρ_dd
                zeros(FT,n_λ,n_layer),      # _ρ_dv
                zeros(FT,n_λ,n_layer),      # _τ_dd
                zeros(FT,n_λ,n_layer),      # _τ_dv
                0                           # _τ_vv
    )
);
