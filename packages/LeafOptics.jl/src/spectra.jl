#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Oct-01: rename the function to leaf_spectra! as the function updates not only fluorescence but also reflectance, transmittance, and absorption spectra
#     2021-Oct-22: add another method to prescribe leaf spectra such as transmittance and reflectance from broadband method
#
#######################################################################################################################################################################################################
"""
This function updates leaf level reflectance, transmittance, and fluorescence spectra related parameters. Supported methods are
- Update leaf spectra based on pigment concentrations
- Update leaf spectra (reflectance and transmittance) to given broadband values
- Update leaf spectra based on pigment concentrations for the entire SPAC

"""
function leaf_spectra! end


#######################################################################################################################################################################################################
#
# Changes made to this method
# General
#     2020-Mar-30: account for carotenoid absorption as PPAR as well as chlorophyll
#     2020-Mar-31: use 40° rather than 59° for _τ_α calculation (following PROSPECT-D)
#     2021-Aug-07: replace function `expint` with that from SpecialFunctions
#     2021-Oct-21: add α to input parameters so that one can roll back to 50° for _τ_α calculation
#     2021-Nov-29: separate HyperspectralAbsorption as constant struct
#     2022-Jan-13: use LeafBiophysics directly in the function rather than Leaf
#     2022-Jun-15: rename LeafBiophysics to HyperspectralLeafBiophysics to be more descriptive
#     2022-Jul-22: add lwc to function variable list
#     2022-Jul-28: run leaf_spectra! only if lwc differs from _v_storage
#     2022-Aug-17: add option reabsorb to function to enable/disable SIF reabsorption
# Bug fix
#     2021-Aug-06: If bio.CBC and bio.PRO are not zero, they are accounted for twice in bio.LMA, thus the spectrum from LMA need to subtract the contribution from CBC and PRO
# To do
#     TODO: add References for this methods
#
#######################################################################################################################################################################################################
"""

    leaf_spectra!(
                bio::HyperspectralLeafBiophysics{FT},
                wls::WaveLengthSet{FT},
                lha::HyperspectralAbsorption{FT},
                lwc::FT;
                APAR_car::Bool = true,
                reabsorb::Bool = true,
                α::FT = FT(40)
    ) where {FT<:AbstractFloat}

Update leaf reflectance and transmittance spectra, and fluorescence spectrum matrices, given
- `bio` `HyperspectralLeafBiophysics` type struct that contains leaf biophysical parameters
- `wls` `WaveLengthSet` type struct that contain wave length bins
- `lha` `HyperspectralAbsorption` type struct that contains absorption characteristic curves
- `lwc` Leaf water content `[mol m⁻²]`
- `APAR_car` If true, carotenoid absorption is accounted for in PPAR, default is `true`
- `reabsorb` If true, SIF reabsorption is enabled; otherwise, mat_b and mat_f should be based on the case with no reabsorption
- `α` Optimum angle of incidence (default is 40° as in PROSPECT-D, SCOPE uses 59°)

# Examples
```julia
wls = EmeraldNamespace.WaveLengthSet{Float64}();
bio = EmeraldNamespace.HyperspectralLeafBiophysics{Float64}();
lha = EmeraldNamespace.HyperspectralAbsorption{Float64}();
leaf_spectra!(bio, wls, lha, 50.0);
leaf_spectra!(bio, wls, lha, 50.0; APAR_car=false);
leaf_spectra!(bio, wls, lha, 50.0; APAR_car=false, α=59.0);
```

"""
leaf_spectra!(
            bio::HyperspectralLeafBiophysics{FT},
            wls::WaveLengthSet{FT},
            lha::HyperspectralAbsorption{FT},
            lwc::FT;
            APAR_car::Bool = true,
            reabsorb::Bool = true,
            α::FT = FT(40)
) where {FT<:AbstractFloat} = (
    # if leaf water content is the same as the historical value, do nothing
    if lwc == bio._v_storage
        return nothing
    end;

    (; MESOPHYLL_N, NDUB) = bio;
    (; K_ANT, K_BROWN, K_CAB, K_CAR_V, K_CAR_Z, K_CBC, K_H₂O, K_LMA, K_PRO, K_PS, NR) = lha;
    (; IΛ_SIF, IΛ_SIFE, Λ_SIF, Λ_SIFE) = wls;

    # calculate the average absorption feature and relative Cab and Car partitions
    bio.k_all    .= (K_CAB   .* bio.cab .+                          # chlorophyll absorption
                     K_CAR_V .* bio.car .* (1 - bio.f_zeax) .+      # violaxanthin carotenoid absorption
                     K_CAR_Z .* bio.car .* bio.f_zeax .+            # zeaxanthin carotenoid absorption
                     K_ANT   .* bio.ant .+                          # anthocynanin absorption absorption
                     K_BROWN .* bio.brown .+                        # TODO: needs to be a concentration
                     K_H₂O   .* (lwc * M_H₂O() / ρ_H₂O() * 100) .+  # water absorption
                     K_CBC   .* bio.cbc .+                          # carbon-based constituents absorption
                     K_PRO   .* bio.pro .+                          # protein absorption
                     K_LMA   .* (bio.lma - bio.cbc - bio.pro)       # dry mass absorption (if some remained)
                    ) ./ MESOPHYLL_N;
    bio.α_cab    .= (K_CAB .* bio.cab) ./ bio.k_all ./ MESOPHYLL_N;
    bio.α_cabcar .= (K_CAB .* bio.cab .+ K_CAR_V .* bio.car .* (1 - bio.f_zeax) .+ K_CAR_Z .* bio.car .* bio.f_zeax) ./ bio.k_all ./ MESOPHYLL_N;

    # calculate the reflectance and transmittance at the interfaces of one layer
    bio._τ   .= (1 .- bio.k_all) .* exp.(-bio.k_all) .+ bio.k_all .^ 2 .* expint.(bio.k_all .+ eps(FT));
    bio._τ_α .= average_transmittance.(α, NR);
    bio._ρ_α .= 1 .- bio._τ_α;
    bio._τ₁₂ .= average_transmittance.(FT(90), NR);
    bio._ρ₁₂ .= 1 .- bio._τ₁₂;
    bio._τ₂₁ .= bio._τ₁₂ ./ (NR .^ 2);
    bio._ρ₂₁ .= 1 .- bio._τ₂₁;

    # top surface side
    bio._denom .= 1 .- (bio._τ .* bio._ρ₂₁) .^ 2;
    bio._τ_top .= bio._τ_α .* bio._τ .* bio._τ₂₁ ./ bio._denom;
    bio._ρ_top .= bio._ρ_α .+ bio._ρ₂₁ .* bio._τ .* bio._τ_top;

    # bottom surface side
    bio._τ_bottom .= bio._τ₁₂ .* bio._τ .* bio._τ₂₁ ./ bio._denom;
    bio._ρ_bottom .= bio._ρ₁₂ .+ bio._ρ₂₁ .* bio._τ .* bio._τ_bottom;

    # calculate the reflectance and transmittance at the interfaces of N layer
    bio._d     .= sqrt.((1 .+ bio._ρ_bottom .+ bio._τ_bottom) .* (1 .+ bio._ρ_bottom .- bio._τ_bottom) .* (1 .- bio._ρ_bottom .+ bio._τ_bottom) .* (1 .- bio._ρ_bottom .- bio._τ_bottom));
    bio._ρ²    .= bio._ρ_bottom .^ 2;
    bio._τ²    .= bio._τ_bottom .^ 2;
    bio._a     .= (1 .+ bio._ρ² .- bio._τ² .+ bio._d) ./ (2 .* bio._ρ_bottom);
    bio._b     .= (1 .- bio._ρ² .+ bio._τ² .+ bio._d) ./ (2 .* bio._τ_bottom);
    bio._bⁿ⁻¹  .= bio._b .^ (MESOPHYLL_N - 1);
    bio._b²ⁿ⁻² .= bio._bⁿ⁻¹ .^ 2;
    bio._a²    .= bio._a .^ 2;
    bio._denom .= bio._a² .* bio._b²ⁿ⁻² .- 1;
    bio._ρ_sub .= bio._a .* (bio._b²ⁿ⁻² .- 1) ./ bio._denom;
    bio._τ_sub .= bio._bⁿ⁻¹ .* (bio._a² .- 1) ./ bio._denom;

    # avoid case of zero absorption
    for _i in eachindex(bio._ρ_bottom)
        if bio._ρ_bottom[_i] + bio._τ_bottom[_i] >= 1
            bio._τ_sub[_i] = bio._τ_bottom[_i] / (bio._τ_bottom[_i] + (1 - bio._τ_bottom[_i]) * (MESOPHYLL_N - 1));
            bio._ρ_sub[_i] = 1 - bio._τ_sub[_i];
        end;
    end;

    # reflectance & transmittance of the leaf: combine top layer with next N-1 layers
    bio._denom .= 1 .- bio._ρ_sub .* bio._ρ_bottom;
    bio.τ_sw   .= bio._τ_top .* bio._τ_sub ./ bio._denom;
    bio.ρ_sw   .= bio._ρ_top .+ bio._τ_top .* bio._ρ_sub .* bio._τ_bottom ./ bio._denom;
    bio.α_sw   .= 1 .- bio.τ_sw .- bio.ρ_sw;

    # Doubling method used to calculate fluoresence is now only applied to the part of the leaf where absorption takes place, that is, the part exclusive of the leaf-air interfaces.
    # The reflectance (rho) and transmittance (tau) of this part of the leaf are now determined by "subtracting" the interfaces.
    # CF Note: All of the below takes about 10 times more time than the RT above. Need to rething speed and accuracy. (10nm is bringing it down a lot!)
    bio._ρ_b .= (bio.ρ_sw .- bio._ρ_α) ./ (bio._τ_α .* bio._τ₂₁ .+ (bio.ρ_sw - bio._ρ_α) .* bio._ρ₂₁);
    bio._tt1 .= bio._τ_α .* bio._τ₂₁;
    bio._tt2 .= bio.τ_sw .* (1 .- bio._ρ_b .* bio._ρ₂₁);
    bio._z   .= bio._tt2 ./ bio._tt1;
    bio._tt1 .= bio._ρ_b .- bio._ρ₂₁ .* bio._z .^ 2;
    bio._tt2 .= 1 .- (bio._ρ₂₁.* bio._z) .^ 2;
    bio._ρ   .= max.(0, bio._tt1 ./ bio._tt2);
    bio._tt1 .= 1 .- bio._ρ_b .* bio._ρ₂₁;
    bio._τ   .= bio._tt1 ./ bio._tt2 .* bio._z;

    # Derive Kubelka-Munk s and k
    for _i in eachindex(bio._ρ)
        if bio._ρ[_i] + bio._τ[_i] < 1
            bio._a[_i] = (1 + bio._ρ[_i] ^ 2 - bio._τ[_i] ^ 2 + bio._d[_i]) / (2 * bio._ρ[_i]);
            bio._b[_i] = (1 - bio._ρ[_i] ^ 2 + bio._τ[_i] ^ 2 + bio._d[_i]) / (2 * bio._τ[_i]);
            bio._d[_i] = sqrt((1 + bio._ρ[_i] + bio._τ[_i]) * (1 + bio._ρ[_i] - bio._τ[_i]) * (1 - bio._ρ[_i] + bio._τ[_i]) *  (1 - bio._ρ[_i] - bio._τ[_i]));
        else
            bio._a[_i] = 1;
            bio._b[_i] = 1;
        end;
    end;

    bio._s     .= bio._ρ ./ bio._τ;
    bio._k     .= log.(bio._b);
    for _i in eachindex(bio._a)
        if 1 < bio._a[_i] < Inf
            bio._s[_i] = 2 * bio._a[_i] / (bio._a[_i] ^ 2 - 1) * log(bio._b[_i]);
            bio._k[_i] = (bio._a[_i] - 1) / (bio._a[_i] + 1) * log(bio._b[_i]);
        end;
    end;
    bio._k_chl .= (APAR_car ? bio.α_cabcar : bio.α_cab) .* bio._k;

    # indices of WLE and WLF within wlp
    _ϵ = FT(2) ^ -NDUB;
    bio._ρ_e     .= view(bio._s,IΛ_SIFE) .* _ϵ;
    bio._ρ_f     .= view(bio._s,IΛ_SIF) * _ϵ;
    bio._sigmoid .= 1 ./ (1 .+ exp.(-Λ_SIF ./ 10) .* exp.(Λ_SIFE' ./ 10));
    bio._mat_f   .= K_PS[IΛ_SIF] .* _ϵ ./ 2 .* bio._k_chl[IΛ_SIFE]' .* bio._sigmoid;
    bio._mat_b   .= K_PS[IΛ_SIF] .* _ϵ ./ 2 .* bio._k_chl[IΛ_SIFE]' .* bio._sigmoid;
    bio._τ_e     .= 1 .- (view(bio._k,IΛ_SIFE) .+ view(bio._s,IΛ_SIFE)) .* _ϵ;
    if reabsorb
        bio._τ_f .= 1 .- (view(bio._k,IΛ_SIF) .+ view(bio._s,IΛ_SIF)) .* _ϵ;
    else
        bio._τ_f .= 1 .- view(bio._s,IΛ_SIF) .* _ϵ;
    end;

    # Doubling adding routine
    for _ in 1:NDUB
        bio._x_e     .= bio._τ_e ./ (1 .- bio._ρ_e .^ 2);
        bio._x_f     .= bio._τ_f ./ (1 .- bio._ρ_f .^ 2);
        bio._τ_e_n   .= bio._τ_e .* bio._x_e;
        bio._τ_f_n   .= bio._τ_f .* bio._x_f;
        bio._ρ_e_n   .= bio._ρ_e .* (1 .+ bio._τ_e_n);
        bio._ρ_f_n   .= bio._ρ_f .* (1 .+ bio._τ_f_n);
        #bio._a₁₁     .= bio._x_f * bio._1_e .+ bio._1_f * bio._x_e';
        #bio._a₁₂     .= (bio._x_f * bio._x_e') .* (bio._ρ_f * bio._1_e .+ bio._1_f * bio._ρ_e');
        #bio._a₂₁     .= 1 .+ (bio._x_f * bio._x_e') .* (1 .+ bio._ρ_f * bio._ρ_e');
        #bio._a₂₂     .= (bio._x_f .* bio._ρ_f) * bio._1_e .+ bio._1_f * (bio._x_e.*bio._ρ_e)';
        bio._a₁₁     .= bio._x_f .* bio._1_e .+ bio._1_f .* bio._x_e';
        bio._a₁₂     .= (bio._x_f .* bio._x_e') .* (bio._ρ_f .* bio._1_e .+ bio._1_f .* bio._ρ_e');
        bio._a₂₁     .= 1 .+ (bio._x_f * bio._x_e') .* (1 .+ bio._ρ_f * bio._ρ_e');
        bio._z_e     .= bio._x_e .* bio._ρ_e;
        bio._z_f     .= bio._x_f .* bio._ρ_f;
        bio._a₂₂     .= bio._z_f .* bio._1_e .+ bio._1_f .* bio._z_e';
        bio._mat_f_n .= bio._mat_f .* bio._a₁₁ .+ bio._mat_b .* bio._a₁₂;
        bio._mat_b_n .= bio._mat_b .* bio._a₂₁ .+ bio._mat_f .* bio._a₂₂;
        bio._τ_e     .= bio._τ_e_n;
        bio._ρ_e     .= bio._ρ_e_n;
        bio._τ_f     .= bio._τ_f_n;
        bio._ρ_f     .= bio._ρ_f_n;
        bio._mat_f   .= bio._mat_f_n;
        bio._mat_b   .= bio._mat_b_n;
    end;

    # This reduced red SIF quite a bit in backscatter, not sure why.
    bio._ρ_b  .= bio._ρ .+ bio._τ .^ 2 .* bio._ρ₂₁ ./ (1 .- bio._ρ .* bio._ρ₂₁);
    bio._z_e  .= view(bio._τ_α,IΛ_SIFE) ./ (1 .- view(bio._ρ₂₁,IΛ_SIFE) .* view(bio._ρ_b,IΛ_SIFE));
    bio._m_xe .= bio._1_f .* bio._z_e';
    bio._m_xf .= view(bio._τ₂₁,IΛ_SIF) ./ (1 .- view(bio._ρ₂₁,IΛ_SIF) .* view(bio._ρ_b,IΛ_SIF)) .* bio._1_e;
    bio._z_e  .= view(bio._τ,IΛ_SIFE) .* view(bio._ρ₂₁,IΛ_SIFE) ./ (1 .- view(bio._ρ,IΛ_SIFE) .* view(bio._ρ₂₁,IΛ_SIFE));
    bio._m_ye .= bio._1_f .* bio._z_e';
    bio._m_yf .= view(bio._τ,IΛ_SIF) .* view(bio._ρ₂₁,IΛ_SIF) ./ (1 .- view(bio._ρ,IΛ_SIF) .* view(bio._ρ₂₁,IΛ_SIF)) .* bio._1_e;
    bio._ma   .= bio._m_xe .* (1 .+ bio._m_ye .* bio._m_yf) .* bio._m_xf;
    bio._mb   .= bio._m_xe .* (bio._m_ye .+ bio._m_yf) .* bio._m_xf;

    bio.mat_b .= bio._ma .* bio._mat_b .+ bio._mb .* bio._mat_f;
    bio.mat_f .= bio._ma .* bio._mat_f .+ bio._mb .* bio._mat_b;

    # store leaf water content
    bio._v_storage = lwc;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes made to this method
# General
#     2021-Oct-22: add another method to prescribe leaf spectra such as transmittance and reflectance from broadband method
#     2021-Nov-29: separate HyperspectralAbsorption as constant struct
#     2022-Jan-13: use LeafBiophysics directly in the function rather than Leaf
#     2022-Jun-15: rename LeafBiophysics to HyperspectralLeafBiophysics to be more descriptive
#
#######################################################################################################################################################################################################
"""

    leaf_spectra!(bio::HyperspectralLeafBiophysics{FT}, wls::WaveLengthSet{FT}, ρ_par::FT, ρ_nir::FT, τ_par::FT, τ_nir::FT) where {FT<:AbstractFloat}

Update leaf reflectance and transmittance (e.g., prescribe broadband PAR and NIR values), given
- `bio` `HyperspectralLeafBiophysics` type struct that contains leaf biophysical parameters
- `wls` `WaveLengthSet` type struct that contain wave length bins
- `ρ_par` Reflectance at PAR region
- `ρ_nir` Reflectance at NIR region
- `τ_par` Transmittance at PAR region
- `τ_nir` Transmittance at NIR region

# Examples
```julia
wls = EmeraldNamespace.WaveLengthSet{Float64}();
bio = EmeraldNamespace.HyperspectralLeafBiophysics{Float64}();
leaf_spectra!(bio, wls, 0.1, 0.45, 0.05, 0.25);
```

"""
leaf_spectra!(bio::HyperspectralLeafBiophysics{FT}, wls::WaveLengthSet{FT}, ρ_par::FT, ρ_nir::FT, τ_par::FT, τ_nir::FT) where {FT<:AbstractFloat} = (
    (; IΛ_NIR, IΛ_PAR) = wls;

    bio.ρ_sw[IΛ_PAR] .= ρ_par;
    bio.ρ_sw[IΛ_NIR] .= ρ_nir;
    bio.τ_sw[IΛ_PAR] .= τ_par;
    bio.τ_sw[IΛ_NIR] .= τ_nir;

    bio.α_sw = 1 .- bio.τ_sw .- bio.ρ_sw;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes made to this method
# General
#     2022-Jun-29: add method for MonoMLGrassSPAC, MonoMLPalmSPAC, and MonoMLTreeSPAC
#
#######################################################################################################################################################################################################
"""

    leaf_spectra!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat}

Update leaf reflectance and transmittance for SPAC, given
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` type SPAC

"""
leaf_spectra!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat} = (
    (; CANOPY, LEAVES) = spac;

    for _leaf in LEAVES
        leaf_spectra!(_leaf.BIO, CANOPY.WLSET, CANOPY.LHA, _leaf.HS.v_storage; APAR_car = _leaf.APAR_CAR);
    end;

    return nothing
);
