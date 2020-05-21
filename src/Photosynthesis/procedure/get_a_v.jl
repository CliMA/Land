"""
    get_a_v(; model::AbstractPhotosynthesisModel, v25::FT, p_i::FT, p_O₂::FT, t_leaf::FT)

Gross photosynthetic rate limited by carboxylation, given
- `model` A PhotosynthesisModelVcJ type parameter sets that stores temperature denpendencies
- `v25` Maximal carboxylation rate at 298.15 K
- `p_i` Leaf internal CO₂ partial pressure
- `p_O₂` O₂ partial pressure
- `t_leaf` Leaf temperature
"""
function get_a_v(;
                 model::AbstractPhotoModelParaSet = PMPSVcJBernacchi{FT}(),
                   v25::FT                        = FT(80.0),
                   p_i::FT                        = FT(30.0),
                  p_O₂::FT                        = FT(21278.25),
                t_leaf::FT                        = FT(298.15)) where {FT}
    vmax   = get_vmax(model.VcT, v25, t_leaf)
    kc     = get_kc(model.KcT, t_leaf)
    ko     = get_ko(model.KoT, t_leaf)
    km     = kc * (1 + p_O₂/ko)
    Γ_star = get_Γ_star(model.ΓsT, t_leaf)
    a_v    = vmax * (p_i-Γ_star) / (p_i+km)
    return a_v
end
