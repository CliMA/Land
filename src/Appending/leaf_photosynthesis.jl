abstract type AbstractPhotosynthesisModel end

abstract type AbstractC3Photosynthesis <: AbstractPhotosynthesisModel end
"""
Classical C3 Photosynthesis model from FvCB, but using a pre-defined gs value as constraint
"""
struct  C3FvCBPhotoGs  <: AbstractC3Photosynthesis end





############  model::C3FvCBPhotoGs #################################
"""
    rubisco_limited_rate!(model::C3FvCBPhotoGs, leaf, met)
Solution when Rubisco activity is limiting and a given gs
Solves quadratic equation 
"""
function rubisco_limited_rate!(model::C3FvCBPhotoGs, leaf, met)
    @unpack  Γstar, Cc, Kc, Ko, o₂, Vcmax, gs, gm, Rd, Γstar = leaf
    @unpack  Ca,ra,g_m_s_to_mol_m2_s = met
    d = Kc*(1+met.ppm_to_Pa*1e6*o₂/Ko)
    a = (ra/g_m_s_to_mol_m2_s + 1.6/gs + 1.0/gm) # = 1/gleaf
    b = -(Ca + d) - (Vcmax-Rd)*a
    c = Vcmax * (Ca - Γstar) - Rd * (Ca+d)
    leaf.Ac = lower_quadratic(a,b,c)
end

"""
    light_limited_rate!(model::C3FvCBPhotoGs,  leaf, met, APAR)
Solution when Electron Transport Rate is limiting and a given gs
C3: RuBP-limited photosynthesis (this is the NADPH requirement stochiometry)
Using curvature θ_j and Jmax
"""
function light_limited_rate!(model::C3FvCBPhotoGs,  leaf, met, APAR)
    @unpack  Γstar, Cc, gs, gm, Rd = leaf
    @unpack  Ca, ra, g_m_s_to_mol_m2_s = met
    electron_transport_rate!(model, leaf, APAR)
    # Check out Bonan's book, this equates diffusion limited and demand limited 
    a = 4*(ra/g_m_s_to_mol_m2_s + 1.6/gs + 1.0/gm) # = 1/gleaf
    b = -(4Ca + 8Γstar) - (leaf.Je - 4Rd)*a/4
    c = leaf.Je * (Ca-Γstar) - Rd * (4Ca + 8Γstar)
    leaf.Aj = lower_quadratic(a,b,c)
end
