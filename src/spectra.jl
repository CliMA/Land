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
# Bug fix
#     2021-Aug-06: If bio.CBC and bio.PRO are not zero, they are accounted for twice in bio.LMA, thus the spectrum from LMA need to subtract the contribution from CBC and PRO
# To do
#     TODO: make brown pigment and absoption curve more general using realistic units
#     TODO: add References for this methods
#     TODO: speed up this function by preallocate memories using a cache structure
#
#######################################################################################################################################################################################################
"""

    leaf_spectra!(bio::HyperspectralLeafBiophysics{FT}, wls::WaveLengthSet{FT}, lha::HyperspectralAbsorption{FT}, lwc::FT; APAR_car::Bool = true, α::FT=FT(40)) where {FT<:AbstractFloat}

Update leaf reflectance and transmittance spectra, and fluorescence spectrum matrices, given
- `bio` `ClimaCache.HyperspectralLeafBiophysics` type struct that contains leaf biophysical parameters
- `wls` `ClimaCache.WaveLengthSet` type struct that contain wave length bins
- `lha` `ClimaCache.HyperspectralAbsorption` type struct that contains absorption characteristic curves
- `lwc` Leaf water content `[mol m⁻²]`
- `APAR_car` If true, carotenoid absorption is accounted for in PPAR, default is `true`
- `α` Optimum angle of incidence (default is 40° as in PROSPECT-D, SCOPE uses 59°)

# Examples
```julia
wls = WaveLengthSet{Float64}();
bio = HyperspectralLeafBiophysics{Float64}(wls);
lha = HyperspectralAbsorption{Float64}(wls);
leaf_spectra!(bio, wls, lha, 50.0);
leaf_spectra!(bio, wls, lha, 50.0; APAR_car=false);
leaf_spectra!(bio, wls, lha, 50.0; APAR_car=false, α=59.0);
```

"""
leaf_spectra!(bio::HyperspectralLeafBiophysics{FT}, wls::WaveLengthSet{FT}, lha::HyperspectralAbsorption{FT}, lwc::FT; APAR_car::Bool = true, α::FT=FT(40)) where {FT<:AbstractFloat} = (
    @unpack MESOPHYLL_N, NDUB = bio;
    @unpack K_ANT, K_BROWN, K_CAB, K_CAR_V, K_CAR_Z, K_CBC, K_H₂O, K_LMA, K_PRO, K_PS, NR = lha;
    @unpack IΛ_SIF, IΛ_SIFE, Λ_SIF, Λ_SIFE = wls;

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
    _τ   = (1 .- bio.k_all) .* exp.(-bio.k_all) .+ bio.k_all .^ 2 .* expint.(bio.k_all .+ eps(FT));
    _τ_α = average_transmittance.(α, NR);
    _ρ_α = 1 .- _τ_α;
    _τ₁₂ = average_transmittance.(FT(90), NR);
    _ρ₁₂ = 1 .- _τ₁₂;
    _τ₂₁ = _τ₁₂ ./ (NR .^ 2);
    _ρ₂₁ = 1 .- _τ₂₁;

    # top surface side
    _denom = 1 .- (_τ .* _ρ₂₁) .^ 2;
    _τ_top = _τ_α .* _τ .* _τ₂₁ ./ _denom;
    _ρ_top = _ρ_α .+ _ρ₂₁ .* _τ .* _τ_top;

    # bottom surface side
    _τ_bottom = _τ₁₂ .* _τ .* _τ₂₁ ./ _denom;
    _ρ_bottom = _ρ₁₂ .+ _ρ₂₁ .* _τ .* _τ_bottom;

    # calculate the reflectance and transmittance at the interfaces of N layer
    _d     = sqrt.((1 .+ _ρ_bottom .+ _τ_bottom) .* (1 .+ _ρ_bottom .- _τ_bottom) .* (1 .- _ρ_bottom .+ _τ_bottom) .* (1 .- _ρ_bottom .- _τ_bottom));
    _ρ²    = _ρ_bottom .^ 2;
    _τ²    = _τ_bottom .^ 2;
    _a     = (1 .+ _ρ² .- _τ² .+ _d) ./ (2 .* _ρ_bottom);
    _b     = (1 .- _ρ² .+ _τ² .+ _d) ./ (2 .* _τ_bottom);
    _bⁿ⁻¹  = _b .^ (MESOPHYLL_N - 1);
    _b²ⁿ⁻² = _bⁿ⁻¹ .^ 2;
    _a²    = _a .^ 2;
    _denom = _a² .* _b²ⁿ⁻² .- 1;
    _ρ_sub = _a .* (_b²ⁿ⁻² .- 1) ./ _denom;
    _τ_sub = _bⁿ⁻¹ .* (_a² .- 1) ./ _denom;

    # avoid case of zero absorption
    _j = findall(_ρ_bottom .+ _τ_bottom .>= 1);
    _τ_sub[_j] = _τ_bottom[_j] ./ (_τ_bottom[_j] + (1 .- _τ_bottom[_j]) * (MESOPHYLL_N - 1));
    _ρ_sub[_j] = 1 .- _τ_sub[_j];

    # reflectance & transmittance of the leaf: combine top layer with next N-1 layers
    _denom   = 1 .- _ρ_sub .* _ρ_bottom;
    bio.τ_sw = _τ_top .* _τ_sub ./ _denom;
    bio.ρ_sw = _ρ_top .+ _τ_top .* _ρ_sub .* _τ_bottom ./ _denom;
    bio.α_sw = 1 .- bio.τ_sw .- bio.ρ_sw;

    # Doubling method used to calculate fluoresence is now only applied to the part of the leaf where absorption takes place, that is, the part exclusive of the leaf-air interfaces.
    # The reflectance (rho) and transmittance (tau) of this part of the leaf are now determined by "subtracting" the interfaces.
    # CF Note: All of the below takes about 10 times more time than the RT above. Need to rething speed and accuracy. (10nm is bringing it down a lot!)
    _ρ_b = (bio.ρ_sw .- _ρ_α) ./ (_τ_α .* _τ₂₁ .+ (bio.ρ_sw - _ρ_α) .* _ρ₂₁);
    _tt1 = _τ_α .* _τ₂₁;
    _tt2 = bio.τ_sw .* (1 .- _ρ_b .* _ρ₂₁);
    _z   = _tt2 ./ _tt1;
    _tt1 = _ρ_b - _ρ₂₁ .* _z .^ 2;
    _tt2 = 1 .- (_ρ₂₁.* _z) .^ 2;
    _ρ   = max.(0, _tt1 ./ _tt2);
    _tt1 = 1 .- _ρ_b .* _ρ₂₁;
    _τ   = _tt1 ./ _tt2 .* _z;

    # Derive Kubelka-Munk s and k
    _i      = findall((_ρ .+ _τ) .< 1);
    _j      = findall((_ρ .+ _τ) .> 1);
    _d[_i] .= sqrt.((1 .+ _ρ[_i] .+ _τ[_i]) .* (1 .+ _ρ[_i] .- _τ[_i]) .* (1 .- _ρ[_i] .+ _τ[_i]) .*  (1 .- _ρ[_i] .- _τ[_i]));
    _a[_i] .= (1 .+ _ρ[_i] .^ 2 .- _τ[_i] .^ 2 .+ _d[_i]) ./ (2 .* _ρ[_i]);
    _b[_i] .= (1 .- _ρ[_i] .^ 2 .+ _τ[_i] .^ 2 .+ _d[_i]) ./ (2 .* _τ[_i]);
    _a[_j] .= 1;
    _b[_j] .= 1;

    _i      = findall((_a .> 1) .& (_a .!= Inf));
    _s      = _ρ ./ _τ;
    _k      = log.(_b);
    _s[_i] .= 2 .* _a[_i] ./ (_a[_i] .^ 2 .- 1) .* log.(_b[_i]);
    _k[_i] .= (_a[_i] .- 1) ./ (_a[_i] .+ 1) .* log.(_b[_i]);
    _k_chl  = (APAR_car ? bio.α_cabcar : bio.α_cab) .* _k;

    # indices of WLE and WLF within wlp
    _ϵ       = FT(2) ^ -NDUB;
    _τ_e     = 1 .- (_k[IΛ_SIFE] .+ _s[IΛ_SIFE]) * _ϵ;
    _τ_f     = 1 .- (_k[IΛ_SIF] .+ _s[IΛ_SIF]) * _ϵ;
    _ρ_e     = _s[IΛ_SIFE] * _ϵ;
    _ρ_f     = _s[IΛ_SIF] * _ϵ;
    _sigmoid = 1 ./ (1 .+ exp.(-Λ_SIF ./ 10) .* exp.(Λ_SIFE' ./ 10));
    _mat_f   = K_PS[IΛ_SIF] .* _ϵ ./ 2 .* _k_chl[IΛ_SIFE]' .* _sigmoid;
    _mat_b   = K_PS[IΛ_SIF] .* _ϵ ./ 2 .* _k_chl[IΛ_SIFE]' .* _sigmoid;

    # Doubling adding routine
    _1_h = ones(FT, 1, length(_τ_e));
    _1_v = ones(FT, length(_τ_f), 1);
    for i in 1:NDUB
        _x_e     = _τ_e ./ (1 .- _ρ_e .^ 2);
        _x_f     = _τ_f ./ (1 .- _ρ_f .^ 2);
        _τ_e_n   = _τ_e .* _x_e;
        _τ_f_n   = _τ_f .* _x_f;
        _ρ_e_n   = _ρ_e .* (1 .+ _τ_e_n);
        _ρ_f_n   = _ρ_f .* (1 .+ _τ_f_n);
        _a₁₁     = _x_f * _1_h .+ _1_v * _x_e';
        _a₁₂     = (_x_f * _x_e') .* (_ρ_f * _1_h .+ _1_v * _ρ_e');
        _a₂₁     = 1 .+ (_x_f * _x_e') .* (1 .+ _ρ_f * _ρ_e');
        _a₂₂     = (_x_f .* _ρ_f) * _1_h .+ _1_v * (_x_e.*_ρ_e)';
        _mat_f_n = _mat_f .* _a₁₁ .+ _mat_b .* _a₁₂;
        _mat_b_n = _mat_b .* _a₂₁ .+ _mat_f .* _a₂₂;
        _τ_e     = _τ_e_n;
        _ρ_e     = _ρ_e_n;
        _τ_f     = _τ_f_n;
        _ρ_f     = _ρ_f_n;
        _mat_f   = _mat_f_n;
        _mat_b   = _mat_b_n;
    end;

    # This reduced red SIF quite a bit in backscatter, not sure why.
    _ρ_b = _ρ .+ _τ .^ 2 .* _ρ₂₁ ./ (1 .- _ρ .* _ρ₂₁);
    _x_e = _1_v * (_τ_α[IΛ_SIFE] ./ (1 .- _ρ₂₁[IΛ_SIFE] .* _ρ_b[IΛ_SIFE]))';
    _x_f = _τ₂₁[IΛ_SIF] ./ (1 .- _ρ₂₁[IΛ_SIF] .* _ρ_b[IΛ_SIF]) * _1_h;
    _y_e = _1_v * (_τ[IΛ_SIFE] .* _ρ₂₁[IΛ_SIFE] ./ (1 .- _ρ[IΛ_SIFE] .* _ρ₂₁[IΛ_SIFE]))';
    _y_f = _τ[IΛ_SIF] .* _ρ₂₁[IΛ_SIF] ./ (1 .- _ρ[IΛ_SIF] .* _ρ₂₁[IΛ_SIF]) * _1_h;
    _a   = _x_e .* (1 .+ _y_e .* _y_f) .* _x_f;
    _b   = _x_e .* (_y_e .+ _y_f) .* _x_f;

    bio.mat_b = _a .* _mat_b + _b .* _mat_f;
    bio.mat_f = _a .* _mat_f + _b .* _mat_b;

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
- `bio` `ClimaCache.HyperspectralLeafBiophysics` type struct that contains leaf biophysical parameters
- `wls` `ClimaCache.WaveLengthSet` type struct that contain wave length bins
- `ρ_par` Reflectance at PAR region
- `ρ_nir` Reflectance at NIR region
- `τ_par` Transmittance at PAR region
- `τ_nir` Transmittance at NIR region

# Examples
```julia
wls = WaveLengthSet{Float64}();
bio = HyperspectralLeafBiophysics{Float64}(wls);
leaf_spectra!(bio, wls, 0.1, 0.45, 0.05, 0.25);
```

"""
leaf_spectra!(bio::HyperspectralLeafBiophysics{FT}, wls::WaveLengthSet{FT}, ρ_par::FT, ρ_nir::FT, τ_par::FT, τ_nir::FT) where {FT<:AbstractFloat} = (
    @unpack IΛ_NIR, IΛ_PAR = wls;

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
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, and `MonoMLTreeSPAC` type SPAC

"""
leaf_spectra!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat} = (
    @unpack LEAVES, LHA, WLSET = spac;

    for _leaf in LEAVES
        leaf_spectra!(_leaf.BIO, WLSET, LHA, _leaf.HS.v_storage; APAR_car = _leaf.APAR_CAR);
    end;

    return nothing
);
