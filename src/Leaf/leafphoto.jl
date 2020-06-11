###############################################################################
#
# Update the temperature dependencies in the leaf
# These functions passed the FT test
# These functions are documented in the Leaf page
#
###############################################################################
"""
    update_leaf_TD!(paraset::C3ParaSet, leaf::LeafParams{FT})
    update_leaf_TD!(paraset::C4ParaSet, leaf::LeafParams{FT})

Update leaf photosynthesis-related temperature dependencies, given
- `paraset` A `C3ParaSet` or `C4ParaSet` type parameter set
- `leaf` A `LeafParams` struct
"""
function update_leaf_TD!(paraset::C3ParaSet, leaf::LeafParams{FT}) where {FT}
    @unpack JT, ReT, VcT, ΓsT = paraset
    @unpack Jmax25, Rd25, T, Vcmax25 = leaf
    KcT = paraset.KcT
    KoT = paraset.KoT
    KpT = paraset.KoT
    ΓsT = paraset.ΓsT

    # update the leaf information
    leaf.Jmax  = Jmax25     * arrhenius_correction( JT, T)
    leaf.Kc    = KcT.VAL_25 * arrhenius_correction(KcT, T)
    leaf.Ko    = KoT.VAL_25 * arrhenius_correction(KoT, T)
    leaf.Rd    = Rd25       * arrhenius_correction(ReT, T)
    leaf.Vcmax = Vcmax25    * arrhenius_correction(VcT, T)
    leaf.Γstar = ΓsT.VAL_25 * arrhenius_correction(ΓsT, T)
end

function update_leaf_TD!(paraset::C4ParaSet, leaf::LeafParams{FT}) where {FT}
    @unpack ReT, VcT, VpT = paraset
    @unpack Rd25, T, Vcmax25, Vpmax25 = leaf
    KpT = paraset.KpT

    # update the leaf information
    leaf.Kpep  = KpT.VAL_25 * arrhenius_correction(KpT, T)
    leaf.Rd    = Rd25       * arrhenius_correction(ReT, T)
    leaf.Vcmax = Vcmax25    * arrhenius_correction(VcT, T)
    leaf.Vpmax = Vpmax25    * arrhenius_correction(VpT, T)
end








###############################################################################
#
# Update the electron transport in the leaf
# This function passed the FT test
# This function is documented in the Leaf page
#
###############################################################################

# TODO add the function for C4 photosynthesis
# TODO if APAR is updated in the leaf, remove the APAR here
"""
    electron_transport_rate!(leaf::LeafParams{FT}, APAR::FT) 

Update the electron transport variables in the leaf struct, given
- `paraset` A `C3ParaSet` type parameter set
- `leaf` A `LeafParams` struct
- `APAR` Absorbed PAR
"""
function electron_transport_rate!(paraset::C3ParaSet, leaf::LeafParams{FT}, APAR::FT) where {FT}
    @unpack maxPSII, Jmax, θ_j, PSII_frac = leaf

    # do the calculations
    _Jp = PSII_frac * maxPSII * APAR
    _b  = _Jp + Jmax
    _c  = _Jp * Jmax
    _J  = ( _b - sqrt(_b^2 - 4*θ_j*_c) ) / (2*θ_j)

    # update the leaf information
    leaf.Je_pot = _Jp
    leaf.Je     = _J
end








###############################################################################
#
# Update the RubisCO-limited photosynthesis in the leaf
# These functions passed the FT test
# These functions are documented in the Leaf page
#
###############################################################################
"""
    rubisco_limited_rate!(paraset::C3ParaSet, leaf::LeafParams{FT})
    rubisco_limited_rate!(paraset::C4ParaSet, leaf::LeafParams{FT})

Update the Rubisco-limited photosynthesis in the leaf struct, given
- `paraset` A `C3ParaSet` or `C4ParaSet` type parameter set
- `leaf` A `LeafParams` struct
"""
function rubisco_limited_rate!(paraset::C3ParaSet, leaf::LeafParams{FT}) where {FT}
    @unpack Γstar, p_i, Kc, Ko, p_O₂, Vcmax = leaf
    leaf.Ac = Vcmax * (p_i - Γstar) / (p_i + Kc*(1 + p_O₂/Ko))
end

function rubisco_limited_rate!(paraset::C4ParaSet, leaf::LeafParams{FT}) where {FT}
    leaf.Ac = leaf.Vcmax
end








###############################################################################
#
# Update the light-limited photosynthesis in the leaf
# These functions passed the FT test
# These functions are documented in the Leaf page
#
###############################################################################
"""
    light_limited_rate!(paraset::C3ParaSet, leaf::LeafParams{FT}, APAR::FT)
    light_limited_rate!(paraset::C4ParaSet, leaf::LeafParams{FT}, APAR::FT)

Update the Electron Transport-limited photosynthesis in the leaf struct, given
- `paraset` A `C3ParaSet` or `C4ParaSet` type parameter set
- `leaf` A `LeafParams` struct
- `APAR` Absorbed PAR
"""
function light_limited_rate!(paraset::C3ParaSet, leaf::LeafParams{FT}, APAR::FT) where {FT}
    # update leaf electron transport
    # TODO This part is indepent of the current function, remove it if Je and Je_pot are updated a priori
    electron_transport_rate!(paraset, leaf, APAR)

    @unpack Γstar, p_i, Je = leaf
    @unpack Eff_1, Eff_2   = paraset

    # do the calculations
    _eff = (p_i - Γstar) / (Eff_1 *p_i + Eff_2*Γstar)
    _Aj  = Je * _eff

    # update the leaf information
    leaf.CO₂_per_electron = _eff
    leaf.Aj               = _Aj
end

function light_limited_rate!(paraset::C4ParaSet, leaf::LeafParams{FT}, APAR::FT) where {FT}
    # do the calculations
    _eff = FT(1/6)
    _Aj  = _eff * (leaf.maxPSII * APAR) / 2

    # update the leaf information
    leaf.CO₂_per_electron = _eff
    leaf.Aj               = _Aj
end








###############################################################################
#
# Update the product-limited photosynthesis in the leaf
# These functions passed the FT test
# These functions are documented in the Leaf page
#
###############################################################################
"""
    product_limited_rate!(paraset::C3ParaSet, leaf::LeafParams{FT})
    product_limited_rate!(paraset::C4ParaSet, leaf::LeafParams{FT})

Update the Product-limited photosynthesis in the leaf struct, given
- `paraset` A `C3ParaSet` or `C4ParaSet` type parameter set
- `leaf` A `LeafParams` struct
"""
function product_limited_rate!(paraset::C3ParaSet, leaf::LeafParams{FT}) where {FT}
    leaf.Ap = leaf.Vcmax / 2
end

function product_limited_rate!(paraset::C4ParaSet, leaf::LeafParams{FT}) where {FT}
    @unpack Kpep, Vpmax, p_i = leaf
    leaf.Ap = max(0, Vpmax * p_i / (p_i + Kpep))
end








###############################################################################
#
# Update the fluorescence in the leaf
# This function passed the FT test
# This function is documented in the Leaf page
#
###############################################################################
"""
    leaf_fluorescence!(paraset::FluoParaSet, leaf::LeafParams{FT})

Compute fluorescence yield, Kn, and Kp for leaf, given
- `paraset` A `FluoParaSet` type parameter set
- `leaf` A `LeafParams` struct
"""
function leaf_fluorescence!(paraset::FluoParaSet, leaf::LeafParams{FT}) where {FT}
    @unpack Kf, Kd, p_i, Γstar, effcon, Ag, maxPSII = leaf
    @unpack Kn1, Kn2, Kn3 = paraset

    #@show leaf.CO2_per_electron
    # Actual effective ETR:
    leaf.Ja = max(0, Ag / leaf.CO₂_per_electron)
    #leaf.Ja = min(leaf.Ja,leaf.Je_pot )
    #@show leaf.Ja

    # Effective photochemical yield:
    if leaf.Ja<= 0
        leaf.φ = maxPSII
    else
        leaf.φ = maxPSII*leaf.Ja/leaf.Je_pot
    end

    #println(flux.Ja, " ", flux.Je_pot)
    leaf.φ   = min(1/maxPSII,leaf.φ)
    # degree of light saturation: 'x' (van der Tol e.Ap. 2014)
    x        = max(0,  1-leaf.φ/leaf.maxPSII)
    #@show x
    #@show leaf.φ
    # Max PSII rate constant
    Kp_max   = FT(4.0)

    x_alpha  = exp(log(x)*Kn2)
    #println(x_alpha)
    
    leaf.Kn  = Kn1 * (1+Kn3)* x_alpha/(Kn3 + x_alpha)
    leaf.Kp  = max(0,-leaf.φ*(Kf+Kd+leaf.Kn)/(leaf.φ-1))

    leaf.Fo  = Kf/(Kf+Kp_max+Kd        )
    leaf.Fo′ = Kf/(Kf+Kp_max+Kd+leaf.Kn)
    leaf.Fm  = Kf/(Kf       +Kd        )
    leaf.Fm′ = Kf/(Kf       +Kd+leaf.Kn)
    leaf.ϕs  = leaf.Fm′*(1-leaf.φ)
    # leaf.eta  = leaf.ϕs/leaf.Fo
    # don't need this anymore, better to use ϕs directly for SIF as Fo is not always fqe=0.01.
    leaf.qQ  = 1-(leaf.ϕs-leaf.Fo′)/(leaf.Fm-leaf.Fo′)
    leaf.qE  = 1-(leaf.Fm-leaf.Fo′)/(leaf.Fm′-leaf.Fo)
    leaf.NPQ = leaf.Kn/(Kf+Kd)
end
