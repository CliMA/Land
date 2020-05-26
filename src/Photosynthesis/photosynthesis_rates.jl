abstract type AbstractPhotosynthesisModel end

abstract type AbstractC3Photosynthesis <: AbstractPhotosynthesisModel end
abstract type AbstractC4Photosynthesis <: AbstractPhotosynthesisModel end

"""
Classical C3 Photosynthesis model from FvCB, using NADPH based stochiometry
"""
struct  C3FvCBPhoto    <: AbstractC3Photosynthesis end
"""
Classical C3 Photosynthesis model from FvCB, but using a pre-defined gs value as constraint
"""
struct  C3FvCBPhotoGs  <: AbstractC3Photosynthesis end
"""
Classical C3 Photosynthesis model from FvCB, using ATP based stochiometry
"""
struct  C3FvCBPhotoATP <: AbstractC3Photosynthesis end
"""
Classical C4 Photosynthesis model from Collatz (but using MM kinetics for PEP-carboxylase)
"""
struct  C4CollatzPhoto <: AbstractC4Photosynthesis end

"""
Compute ETR (Je) from Jmax and Jpot using quadratic
"""
function electron_transport_rate!(model::AbstractPhotosynthesisModel, leaf, APAR)
    @unpack  maxPSII, Jmax,θ_j, PSII_frac = leaf
    leaf.Je_pot = PSII_frac * maxPSII .* APAR;
    lower_quadratic!(θ_j, -(leaf.Je_pot + Jmax), leaf.Je_pot * Jmax, leaf.Je)
end

function electron_transport_rate2!(model::AbstractPhotosynthesisModel, leaf, APAR)
    @unpack  maxPSII, Jmax,θ_j, PSII_frac = leaf
    Je_pot = PSII_frac * maxPSII .* APAR;
    lower_quadratic2!(θ_j, -(Je_pot .+ Jmax), leaf.Je_pot * Jmax, Je_pot)
end

############  model::C3FvCBPhoto #################################
"""
    rubisco_limited_rate(f::C3FvCBPhoto, leaf)
Solution when Rubisco activity is limiting (Classical Farquhar, von Caemmerer, Berry model) 
"""
function rubisco_limited_rate!(model::AbstractC3Photosynthesis, leaf, met)
    @unpack  Γstar, Cc, Kc, Ko, o₂, Vcmax = leaf
    leaf.Ac = Vcmax * (met.ppm_to_Pa*Cc-Γstar) / (met.ppm_to_Pa*Cc + Kc*(1 +met.ppm_to_Pa*1e6*o₂/Ko));
end

"""
    product_limited_rate(f::C3FvCBPhoto, leaf)
Solution when Triose Phosphate Export is limiting (CLM implementation)) 
"""
function product_limited_rate!(model::AbstractC3Photosynthesis, leaf, met)
    leaf.Ap = leaf.Vcmax/2;
end

"""
    light_limited_rate(f::C3FvCBPhoto, flux, leaf)
Solution when Electron Transport Rate is limiting (Classical Farquhar, von Caemmerer, Berry model)
C3: RuBP-limited photosynthesis (this is the NADPH requirement stochiometry)
Using curvature θ_j and Jmax
"""
function light_limited_rate!(model::AbstractC3Photosynthesis, leaf, met, APAR)
    @unpack  Γstar, Cc = leaf
    electron_transport_rate!(model, leaf, APAR)
    leaf.CO2_per_electron = (met.ppm_to_Pa*Cc-Γstar) / (4*met.ppm_to_Pa*Cc + 8Γstar) ;
    leaf.Aj = leaf.Je * leaf.CO2_per_electron;
end

############  model::C3FvCBPhotoATP #################################
"""
    light_limited_rate(f::C3FvCBPhotoATP, leaf)
Solution when Electron Transport Rate is limiting (Classical Farquhar, von Caemmerer, Berry model)
C3: RuBP-limited photosynthesis (this is the ATP requirement stochiometry)
Using curvature θ_j and Jmax
"""
function light_limited_rate!(model::C3FvCBPhotoATP, leaf, met, APAR)
    @unpack  Γstar, Cc = leaf
    electron_transport_rate!(model, leaf, APAR)
    leaf.Aj = leaf.Je * (met.ppm_to_Pa*Cc-Γstar) / (4.5*met.ppm_to_Pa*Cc + 10.5*Γstar)
end


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
    lower_quadratic!(a,b,c,leaf.Ac)
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
    lower_quadratic(a,b,c,leaf.Aj)
end

############  model::C4_Collatz_Photo #################################
"""
    rubisco_limited_rate(f::C4_Collatz_Photo, leaf)
Solution when C4 Rubisco activity is limiting (Collatz) 
"""
function rubisco_limited_rate!(model::AbstractC4Photosynthesis, leaf, met)
    leaf.Ac = leaf.Vcmax
end

"""
    light_limited_rate!(model::C4CollatzPhoto, leaf, met, APAR)
Solution when C4 Electron Transport Rate is limiting (Collatz model)
Linear with J
"""
function light_limited_rate!(model::C4CollatzPhoto, leaf, met, APAR)
    FT = typeof(APAR)
    # Efficiency of ETR for C4 plants ()
    α= FT(1/6)
    leaf.CO2_per_electron = α;
    leaf.Aj = α*(leaf.maxPSII * APAR)/2;
end

"""
    product_limited_rate(f::C4_Collatz_Photo, leaf)
Solution for PEP carboxylase-limited rate of carboxylation
"""
function product_limited_rate!(model::C4CollatzPhoto, leaf, met)
    @unpack  Kpep,Vpmax, Cc = leaf
    leaf.Ap = max(0,Vpmax * (met.ppm_to_Pa*Cc) / (met.ppm_to_Pa*Cc + Kpep))
end











