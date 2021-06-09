###############################################################################
#
# Calculate soil albedo based on color class and soil moisture
#
###############################################################################
function soil_albedos(color::Int, swc::FT) where {FT<:AbstractFloat}
    # soil class must be from 1 to 20
    @assert 1 <= color <=20;

    # calculate albedos for PAR and NIR
    _delta   = max(0, FT(0.11) - FT(0.4) * swc);
    _alb_par = max(SOIL_BNDS[color,1], SOIL_BNDS[color,3] + _delta);
    _alb_nir = max(SOIL_BNDS[color,2], SOIL_BNDS[color,4] + _delta);

    return _alb_par,_alb_nir
end



# TODO add multiple method selections
function fit_soil_mat!(soil::SoilOpticals{FT}, wls::WaveLengths{FT}, swc::FT) where {FT<:AbstractFloat}
    # fit the curve only if the values mismatch
    _ref_par,_ref_nir = soil_albedos(soil.color, swc);
    if soil.ρ_PAR !== _ref_par || soil.ρ_NIR !== _ref_nir
        fit_soil_mat!(soil, wls, _ref_par, _ref_nir);
    end

    return nothing
end




# TODO need to check memory allocations
function fit_soil_mat!(soil::SoilOpticals{FT}, wls::WaveLengths{FT}, ref_PAR::FT, ref_NIR::FT) where {FT<:AbstractFloat}
    @unpack iWLF = wls;
    @unpack dry_NIR, dry_PAR, wet_NIR, wet_PAR = soil;

    # update soil PAR and NIR albedo
    soil.ρ_PAR = ref_PAR;
    soil.ρ_NIR = ref_NIR;

    # solve for weights
    soil.SW_vec_2[1] = (ref_PAR * wet_NIR - ref_NIR * wet_PAR) /
                       (wet_NIR * dry_PAR - wet_PAR * dry_NIR);
    soil.SW_vec_2[2] = (ref_PAR * dry_NIR - ref_NIR * dry_PAR) /
                       (dry_NIR * wet_PAR - dry_PAR * wet_NIR);

    # update vectors in soil
    mul!(soil.ρ_SW, soil.SW_mat_2, soil.SW_vec_2);
    soil.ρ_SW_SIF .= view(soil.ρ_SW, iWLF);

    return nothing
end




# TODO need to check memory allocations
function fit_soil_mat_4!(soil::SoilOpticals{FT}, wls::WaveLengths{FT}, ref_PAR::FT, ref_NIR::FT) where {FT<:AbstractFloat}
    @unpack iWLF = wls;
    @unpack dry_NIR, dry_PAR, wet_NIR, wet_PAR = soil;

    # update soil PAR and NIR albedo
    soil.ρ_PAR = ref_PAR;
    soil.ρ_NIR = ref_NIR;
    # TODO weighted mean?
    @inline _fit(x::Vector{FT}) where {FT<:AbstractFloat} = (
        soil.ρ_SW_raw .= max.(0, soil.SW_mat_raw_4 * x);
        _diff = (soil.ρ_SW_raw[ 25] - ref_PAR)^2 +
                (soil.ρ_SW_raw[ 26] - ref_PAR)^2 +
                (soil.ρ_SW_raw[ 27] - ref_PAR)^2 +
                (soil.ρ_SW_raw[ 61] - ref_NIR)^2 +
                (soil.ρ_SW_raw[ 71] - ref_NIR)^2 +
                (soil.ρ_SW_raw[ 81] - ref_NIR)^2 +
                (soil.ρ_SW_raw[ 91] - ref_NIR)^2 +
                (soil.ρ_SW_raw[101] - ref_NIR)^2 +
                (soil.ρ_SW_raw[111] - ref_NIR)^2 +
                (soil.ρ_SW_raw[121] - ref_NIR)^2 +
                (soil.ρ_SW_raw[131] - ref_NIR)^2 +
                (soil.ρ_SW_raw[141] - ref_NIR)^2 +
                (soil.ρ_SW_raw[151] - ref_NIR)^2 +
                (soil.ρ_SW_raw[161] - ref_NIR)^2 +
                (soil.ρ_SW_raw[171] - ref_NIR)^2 +
                (soil.ρ_SW_raw[181] - ref_NIR)^2 +
                (soil.ρ_SW_raw[191] - ref_NIR)^2 +
                (soil.ρ_SW_raw[201] - ref_NIR)^2 +
                (soil.ρ_SW_raw[211] - ref_NIR)^2 +
                (soil.ρ_SW_raw[211] - ref_NIR)^2;
        #_alb_par = mean( view(soil.ρ_SW_raw,1:31) );
        #_alb_nir = mean( view(soil.ρ_SW_raw,31:length(soil.ρ_SW_raw)) );
        #_diff = (_alb_par - ref_PAR)^2 + (_alb_nir - ref_NIR)^2;
        @show x, _diff;
        return -_diff
    );

    _ms = ReduceStepMethodND{FT}(x_mins = FT[-2,-2,-2,-2],
                                 x_maxs = FT[2,2,2,2],
                                 x_inis = FT[0.3742, -0.0761, 0.0254, -0.0983],
                                 Δ_inis = FT[0.1,0.1,0.1,0.1]);
    _tol = SolutionToleranceND{FT}(FT[0.001,0.001,0.001,0.001], 50);
    _sol = find_peak(_fit, _ms, _tol);
    soil.SW_vec_4 .= _sol;

    # update vectors in soil
    mul!(soil.ρ_SW, soil.SW_mat_4, soil.SW_vec_4);
    soil.ρ_SW_SIF .= view(soil.ρ_SW, iWLF);

    return nothing
end
