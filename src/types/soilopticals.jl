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
    # hyperspectral soil albedo
    "Shortwave albedo matrix with `WL` from [`Wavelengths`](@ref)"
    SW_mat::Matrix{FT}
    "Shortwave albedo matrix from 4 bands, wavelengths are 400:10:2500 nm"
    SW_mat_raw::Matrix{FT}
    "Shortwave albedo weight from 4 bands"
    SW_vec::Vector{FT}
    "Soil surface temperature"
    T::FT
    "Longtwave albedo"
    ρ_LW::Vector{FT}
    "Shortwave albedo that matches `WL` from [`Wavelengths`](@ref)"
    ρ_SW::Vector{FT}
    "Shortwave albedo that matches `WLF` from [`Wavelengths`](@ref)"
    ρ_SW_SIF::Vector{FT}
    "Shortwave albedo for NIR"
    ρ_NIR::FT
    "Shortwave albedo for PAR"
    ρ_PAR::FT
    "Shortwave absorption that equals `1 - ρ_SW`"
    ε_SW::Vector{FT}

    # constructor
    function SoilOpticals{FT}(wls::WaveLengths{FT}) where {FT<:AbstractFloat}
        @unpack nWL, nWLF, sWL = wls;

        # read data that has a 10 nm stepping
        _dat = read_csv(SOIL_GSV);
        _raw = Matrix{FT}(_dat[:,2:end-1]);

        # extend the data to 1 nm stepping by interpolating the data
        _wlr = _dat.WL;
        _wle = collect(400:2500.1);
        _ext = Matrix{FT}(undef, length(_wle), 4);
        for _ie in 1:size(_ext,1)-1
            _wl= _wle[_ie];
            _ir = Int(fld(_wl - 400, 10)) + 1;
            _a = (_wl - _wlr[_ir]) * 0.1;
            _ext[_ie,:] .= (1-_a) .* _raw[_ir,:] .+ _a .* _raw[_ir+1,:];
        end
        _ext[end,:] = _raw[end,:];

        # rescale the data to match the steppings of wavelength set
        _res = Matrix{FT}(undef, nWL, 4);

        # fill in the arrays
        for _i_res in 1:nWL
            _wo = findall( (_wle .>= sWL[_i_res]) .& (_wle .< sWL[_i_res+1]) )
            if length(_wo) == 0
                @warn "Some wavelengths out of bounds $(string(sWL[_i_res]))";
            end;
            _res[_i_res,1] = mean( _ext[_wo,1] );
            _res[_i_res,2] = mean( _ext[_wo,2] );
            _res[_i_res,3] = mean( _ext[_wo,3] );
            _res[_i_res,4] = mean( _ext[_wo,4] );
        end

        return new(_res, _raw, ones(FT,4), T_25(FT), FT[0.1],
                   ones(FT,nWL)*FT(0.2), ones(FT,nWLF)*FT(0.2), FT(0.2),
                   FT(0.2), ones(FT,nWL)*FT(0.8))
    end
end
