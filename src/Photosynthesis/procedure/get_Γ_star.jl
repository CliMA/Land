"""
    get_Γ_star(paraset::AbstractΓStar, t_leaf::FT)

Γ_star at leaf temperature, given
- `para_set` One `AbstractΓStar` type that store temperature correction information
- `t_leaf` Leaf temperature
"""
function get_Γ_star(paraset::AbstractΓStarTD, t_leaf::FT) where {FT}
    return paraset.Γ_star * arrhenius_correction(paraset, t_leaf)
end
