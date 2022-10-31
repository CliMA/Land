"""

    use_clm_td!(photo_set::C3ParaSet{FT}, t_opt::FT) where {FT<:AbstractFloat}

Update the TD for Vcmax, given
- `photo_set` `C3ParaSet` type parameter set
- `t_opt` 10 day average temperature in `[K]`

"""
function use_clm_td!(photo_set::C3ParaSet{FT}, t_opt::FT) where {FT<:AbstractFloat}
    vtd = photo_set.VcT;
    vtd.ΔHa_to_RT25 = 72000.0 / GAS_R() / T₂₅();
    vtd.ΔHd_to_R = 200000.0 / GAS_R();
    vtd.ΔSv_to_R = (668.39 - 1.07 * (t_opt - T₀())) / GAS_R();
    vtd.C = 1 + exp( vtd.ΔSv_to_R - vtd.ΔHd_to_R/T₂₅(FT) );

    jtd = photo_set.JT;
    jtd.ΔHa_to_RT25 = 50000.0 / GAS_R() / T₂₅();
    jtd.ΔHd_to_R = 200000.0 / GAS_R();
    jtd.ΔSv_to_R = (659.70 - 0.75 * (t_opt - T₀())) / GAS_R();
    jtd.C = 1 + exp( jtd.ΔSv_to_R - jtd.ΔHd_to_R/T₂₅(FT) );

    return nothing
end
