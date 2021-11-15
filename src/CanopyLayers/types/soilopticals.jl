###############################################################################
#
# Soil optical parameters
#
###############################################################################
"""
    mutable struct SoilOpticals{FT}

A struct of soil optical parameters

# Fields
$(TYPEDFIELDS)
"""
mutable struct SoilOpticals{FT}
    "Soil surface temperature"
    T::FT
    "Soil color class"
    color::Int

    # broad band showrtwave albedo
    "Shortwave albedo for NIR"
    ρ_NIR::FT
    "Shortwave albedo for PAR"
    ρ_PAR::FT

    # hyperspectral shortwave albedo
    "Shortwave albedo that matches `WL` from [`WaveLengths`](@ref)"
    ρ_SW::Vector{FT}
    "Shortwave albedo that matches `WLF` from [`WaveLengths`](@ref)"
    ρ_SW_SIF::Vector{FT}
    "Shortwave absorption that equals `1 - ρ_SW`"
    ε_SW::Vector{FT}

    # hyperspectral soil albedo
    "Shortwave albedo matrix from 4 bands, wavelengths are 400:10:2500 nm"
    SW_mat_raw_4::Matrix{FT}
    "Shortwave albedo matrix from 2 bands, wavelengths are 400:10:2500 nm"
    SW_mat_raw_2::Matrix{FT}
    "Shortwave albedo matrix from 4 bands with `WL` from [`WaveLengths`](@ref)"
    SW_mat_4::Matrix{FT}
    "Shortwave albedo matrix from 2 bands with `WL` from [`WaveLengths`](@ref)"
    SW_mat_2::Matrix{FT}
    "Shortwave albedo weight from 4 bands"
    SW_vec_4::Vector{FT}
    "Shortwave albedo weight from 2 bands"
    SW_vec_2::Vector{FT}
    "Shortwave albedo, wavelengths are 400:10:2500 nm"
    ρ_SW_raw::Vector{FT}

    # hyperspectral longwave albedo
    "Longtwave albedo"
    ρ_LW::Vector{FT}

    # cache that stores mean band values (used to speed up calculations)
    "Mean value for day band 1 in NIR region"
    dry_NIR::FT
    "Mean value for day band 1 in PAR region"
    dry_PAR::FT
    "Mean value for day band 1 in NIR region"
    wet_NIR::FT
    "Mean value for day band 1 in PAR region"
    wet_PAR::FT

    # constructor
    function SoilOpticals{FT}(wls::WaveLengths{FT}) where {FT<:AbstractFloat}
        @unpack nWL, nWLF, sWL = wls;

        # read data that has a 10 nm stepping
        _dat = read_csv(SOIL_GSV);
        _raw_4 = Matrix{FT}(_dat[:,2:end-1]);

        # extend the data to 1 nm stepping by interpolating the data
        _wlr = _dat.WL;
        _wle = collect(400:2500.1);
        _ext_4 = Matrix{FT}(undef, length(_wle), 4);
        for _ie in 1:size(_ext_4,1)-1
            _wl= _wle[_ie];
            _ir = Int(fld(_wl - 400, 10)) + 1;
            _a = (_wl - _wlr[_ir]) * 0.1;
            _ext_4[_ie,:] .= (1-_a) .* _raw_4[_ir,:] .+ _a .* _raw_4[_ir+1,:];
        end
        _ext_4[end,:] = _raw_4[end,:];

        # rescale the data to match the steppings of wavelength set
        _res_4 = Matrix{FT}(undef, nWL, 4);

        # fill in the arrays
        for _i_res in 1:nWL
            _wo = findall( (_wle .>= sWL[_i_res]) .& (_wle .< sWL[_i_res+1]) )
            if length(_wo) == 0
                @warn "Some wavelengths out of bounds $(string(sWL[_i_res]))";
            end;
            _res_4[_i_res,1] = mean( _ext_4[_wo,1] );
            _res_4[_i_res,2] = mean( _ext_4[_wo,2] );
            _res_4[_i_res,3] = mean( _ext_4[_wo,3] );
            _res_4[_i_res,4] = mean( _ext_4[_wo,4] );
        end

        # 2 band values
        _raw_2 = [_raw_4[:,1] _raw_4[:,end]];
        _res_2 = [_res_4[:,1] _res_4[:,end]];

        # mean values for 2 bands
        _dry_nir = mean( _raw_2[31:end,1] );
        _dry_par = mean( _raw_2[1:31,1]   );
        _wet_nir = mean( _raw_2[31:end,2] );
        _wet_par = mean( _raw_2[1:31,2]   );

        return new(T_25(FT), 1, FT(0.2), FT(0.2), ones(FT,nWL)*FT(0.2),
                   ones(FT,nWLF)*FT(0.2), ones(FT,nWL)*FT(0.8), _raw_4, _raw_2,
                   _res_4, _res_2, ones(FT,4), ones(FT,2), ones(FT,211),
                   FT[0.1], _dry_nir, _dry_par, _wet_nir, _wet_par)
    end
end
