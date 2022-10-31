#
#
# This file is meant for CliMA Land example
#
#
using DataFrames: DataFrame, DataFrameRow
using Dates: isleapyear
using JLD2: load

using NetcdfIO: read_nc, save_nc!
using PkgUtility: month_days, nanmean

using Land.CanopyLayers: EVI, FourBandsFittingHybrid, NDVI, NIRv, SIF_WL, SIF_740, fit_soil_mat!
using Land.Photosynthesis: C3CLM, use_clm_td!
using Land.PlantHydraulics: VanGenuchten, create_tree
using Land.SoilPlantAirContinuum: CNPP, GPP, PPAR, SPACMono, T_VEG, initialize_spac_canopy!, prescribe_air!, prescribe_swc!, prescribe_t_leaf!, spac_beta_max, update_Cab!, update_LAI!, update_VJRWW!,
      update_par!, update_sif!, zenith_angle
using Land.StomataModels: BetaGLinearPsoil, ESMMedlyn, GswDrive, gas_exchange!, gsw_control!, prognostic_gsw!


DF_VARIABLES  = ["F_H2O", "F_CO2", "F_GPP", "SIF683", "SIF740", "SIF757", "SIF771", "NDVI", "EVI", "NIRv"];


"""

    prepare_wd(dict::Dict, wd_file::String)

Prepare weather driver dataframe to feed CliMA Land, given
- `dict` Dictionary that store grid information
- `wd_file` Weather driver file

"""
function prepare_wd(dict::Dict, wd_file::String)
    _df_in = read_nc(wd_file);

    # compute T_MEAN based on the weather driver
    _df_in[!,"CO2"]         .= 0.0;
    _df_in[!,"Chlorophyll"] .= 0.0;
    _df_in[!,"LAI"]         .= 0.0;
    _df_in[!,"Vcmax"]       .= 0.0;
    _df_in[!,"T_MEAN"]      .= 0.0;
    for _i in eachindex(_df_in.T_MEAN)
        if _i < 240
            _df_in[_i,"T_MEAN"] = nanmean( max.(_df_in.T_AIR[1:_i], _df_in.T_LEAF[1:_i]) );
        else
            _df_in[_i,"T_MEAN"] = nanmean( max.(_df_in.T_AIR[_i-239:_i], _df_in.T_LEAF[_i-239:_i]) );
        end;
    end;

    #
    # extropolate the time series based on input variable dimensions, dimensions must be within supported settings
    #     1. extropolate the data to 1D resolution
    #     2. extropolate the data to 1H resolution
    #
    _year = dict["year"];
    _days = isleapyear(_year) ? 366 : 365;
    @inline nt_to_1h(label::String) = (
        _dat_in = dict[label];
        @assert length(_dat_in) in [366, 365, 53, 52, 46, 12, 1] "Dataset length not supported";

        if length(_dat_in) == 1
            _dat_1d = repeat([_dat_in;]; inner = _days);
        elseif length(_dat_in) == 12
            _dat_1d = [([repeat(_dat_in[_m:_m], month_days(_year, _m)) for _m in 1:12]...)...]
        elseif length(_dat_in) == 46
            _dat_1d = repeat(_dat_in; inner = 8)[1:_days]
        elseif length(_dat_in) in [52,53]
            _dat_1d = repeat([_dat_in;_dat_in[end]]; inner = 7)[1:_days]
        elseif length(_dat_in) in [365,366]
            _dat_1d = [_dat_in;_dat_in[end]][1:_days]
        end;

        return repeat(_dat_1d; inner = 24)
    );
    _df_in[!,"CO2"]         .= nt_to_1h("co2_concentration");
    _df_in[!,"Chlorophyll"] .= nt_to_1h("chlorophyll");
    _df_in[!,"CI"]          .= nt_to_1h("clumping_index");
    _df_in[!,"LAI"]         .= nt_to_1h("leaf_area_index");
    _df_in[!,"Vcmax"]       .= nt_to_1h("vcmax");

    # add the fields to store outputs
    for _label in DF_VARIABLES
        _df_in[!,_label] .= 0.0;
    end;

    return _df_in
end



"""

    prepare_spac(dict::Dict; FT = Float64)

Create a SPAC, given
- `dict` Dictionary of GriddingMachine data in a grid
- `FT` Floating number type

"""
function prepare_spac(dict::Dict; FT = Float64)
    # read general information from dict
    _lat = dict["latitude"];
    _lon = dict["longitude"]
    _sm = ESMMedlyn{FT}();

    # use JULES soil depth 0.00 -- 0.10 -- 0.35 -- 1.00 -- 3.00 m, and assume 2 m deep root (z_root = -2) for all the sites
    _soil_bounds = FT[0, -0.1, -0.35, -1, -3];
    _z_canopy    = max(FT(0.1), dict["canopy_height"]);
    _Δz          = _z_canopy / 20;
    _air_bounds  = collect(0:_Δz:_z_canopy+2*_Δz);
    _plant_hs    = create_tree(FT(-2), _z_canopy/2, _z_canopy, _soil_bounds, _air_bounds);

    # create a SPACMono struct, redefine the wavelength limits for PAR if ePAR is true
    _node = SPACMono{FT}(soil_bounds=_soil_bounds, air_bounds=_air_bounds, z_canopy=_z_canopy, z_root=-2, plant_hs=_plant_hs, latitude=_lat, longitude=_lon, stomata_model=_sm);

    for _iPS in _node.plant_ps
        _iPS.g_min   = eps(FT);
        _iPS.g_min25 = eps(FT);
        _iPS.g_max   = 0.8;
        _iPS.g_max25 = 0.8;
    end;

    # update soil type information per layer
    for _i in eachindex(_node.plant_hs.roots)
        _α  = dict["soil_vg_α"][_i];
        _n  = dict["soil_vg_n"][_i];
        _Θr = dict["soil_vg_Θr"][_i];
        _Θs = dict["soil_vg_Θs"][_i];
        _node.plant_hs.roots[_i].sh = VanGenuchten{FT}(stype = "JULES", α = _α, n = _n, Θs = _Θs, Θr = _Θr);
    end;

    # update leaf mass per area (from m² kg⁻¹ to g cm⁻²)
    _lma = dict["leaf_mass_per_area"];
    for leaf in _node.leaves_rt
        leaf.Cm = _lma;
    end;

    # set up empirical model
    if typeof(_sm) <: ESMMedlyn
        _node.photo_set = C3CLM(FT);
        _node.stomata_model.g1 = dict["g1_medlyn_c3"];
        _node.stomata_model.g0 = 1e-3;
    else
        @warn "Stomatal model parameters are not initialized for $(typeof(_sm))";
    end;

    # update soil color class from CLM dataset
    _node.soil_opt.color = dict["soil_color"];

    # update the Vcmax, Jmax, and Vpmax
    update_VJRWW!(_node, nanmean(dict["vcmax"]));

    # initialize the canopy RT model
    initialize_spac_canopy!(_node);

    return _node
end


"""

Structure that store memory information

"""
Base.@kwdef mutable struct SPACMemory{FT<:AbstractFloat}
    chl::FT = -9999
    lai::FT = -9999
    vcm::FT = -9999
end


"""

    prescribe_parameters!(spac::SPACMono{FT}, dfr::DataFrame, mem::SPACMemory{FT}, deepcopies::Vector) where {FT<:AbstractFloat}

Prescibe parameters for the SPAC, given
- `spac` Soil plant air continuum struct
- `dfr` Weather driver dataframe row
- `mem` Memory cache struct
- `deepcopies` Deepcopies of radiation used to scale direct and diffuse radiation

"""
function prescribe_parameters!(spac::SPACMono{FT}, dfr::DataFrameRow, mem::SPACMemory{FT}, deepcopies::Vector) where {FT<:AbstractFloat}
    # read the data out of dataframe row to reduce memory allocation
    _df_atm::FT = dfr.P_ATM;
    _df_chl::FT = dfr.Chlorophyll;
    _df_cli::FT = dfr.CI;
    _df_co2::FT = dfr.CO2;
    _df_dif::FT = dfr.RAD_DIF;
    _df_dir::FT = dfr.RAD_DIR;
    _df_doy::FT = dfr.FDOY;
    _df_lai::FT = dfr.LAI;
    _df_sw1::FT = dfr.SWC_1;
    _df_sw2::FT = dfr.SWC_2;
    _df_sw3::FT = dfr.SWC_3;
    _df_sw4::FT = dfr.SWC_4;
    _df_tar::FT = dfr.T_AIR;
    _df_tlf::FT = dfr.T_LEAF;
    _df_tmn::FT = dfr.T_MEAN;
    _df_vcm::FT = dfr.Vcmax;
    _df_vpd::FT = dfr.VPD;
    _df_wnd::FT = dfr.WIND;

    # adjust optimum t based on 10 day moving average skin temperature
    use_clm_td!(spac.photo_set, _df_tmn);

    # if total LAI, Vcmax, or Chl changes, update them (add vertical Vcmax profile as well)
    _trigger_lai::Bool = !isnan(_df_lai) && (_df_lai != mem.lai);
    _trigger_vcm::Bool = !isnan(_df_vcm) && (_df_vcm != mem.vcm);
    _trigger_chl::Bool = !isnan(_df_chl) && (_df_chl != mem.chl);
    if _trigger_lai
        update_LAI!(spac, _df_lai);
        mem.lai = _df_lai;
    end;

    if _trigger_lai || _trigger_vcm
        update_VJRWW!(spac, _df_vcm; expo = FT(0.3));
        mem.vcm = _df_vcm;
    end;

    if _trigger_chl
        update_Cab!(spac, _df_chl; cab_2_car = FT(1/7));
        mem.chl = _df_chl;
    end;

    # update clumping index
    spac.canopy_rt.Ω = _df_cli;
    spac.canopy_rt.clump_a = _df_cli;

    # sync the environmental conditions per layer
    prescribe_air!(spac, _df_co2, _df_atm, _df_tar, _df_vpd, _df_wnd);
    prescribe_t_leaf!(spac, max(_df_tar, _df_tlf));

    # run the chunks below only when total radiation is higher than 10
    if _df_dir + _df_dif < 10
        return nothing
    end;

    # update soil water matrices per layer
    prescribe_swc!(spac, _df_sw1, _df_sw2, _df_sw3, _df_sw4);

    # update soil albedo using FourBandsFittingHybrid
    _method = FourBandsFittingHybrid();
    fit_soil_mat!(spac.soil_opt, spac.wl_set, spac.swc[1], _method);

    # update PAR related information
    spac.in_rad.E_direct  .= deepcopies[1].E_direct  .* _df_dir ./ deepcopies[2];
    spac.in_rad.E_diffuse .= deepcopies[1].E_diffuse .* _df_dif ./ deepcopies[3];
    spac.angles.sza = min(88, zenith_angle(spac.latitude, _df_doy));
    update_par!(spac);

    return nothing
end


"""

    update_gsw!(spac::SPACMono{FT}, sm::ESMMedlyn{FT}, ind::Int, δt::FT; β::FT = FT(1)) where {FT<:AbstractFloat}

Wrapper function to use prognostic_gsw!, given
- `spac` Soil plant air continuum struct
- `sm` Medlyn stomtal model
- `ind` Canopy layer number
- `δt` Time step
- `β` Tuning factor

"""
function update_gsw!(spac::SPACMono{FT}, sm::ESMMedlyn{FT}, ind::Int, δt::FT; β::FT = FT(1)) where {FT<:AbstractFloat}
    prognostic_gsw!(spac.plant_ps[ind], spac.envirs[ind], sm, β, δt);

    return nothing
end


"""

    run_time_step!(spac::SPACMono{FT}, dfr::DataFrame) where {FT<:AbstractFloat}

Run CliMA Land in a time step, given
- `spac` Soil plant air continuum struct
- `dfr` Weather driver dataframe row
- `ind` Time index

"""
function run_time_step!(spac::SPACMono{FT}, dfr::DataFrameRow, beta::BetaGLinearPsoil{FT}) where {FT<:AbstractFloat}
    # read the data out of dataframe row to reduce memory allocation
    _df_dif::FT = dfr.RAD_DIF;
    _df_dir::FT = dfr.RAD_DIR;

    # compute beta factor (based on Psoil, so canopy does not matter)
    _βm = spac_beta_max(spac, beta);

    # calculate leaf level flux per canopy layer
    for _i_can in 1:spac.n_canopy
        _iEN = spac.envirs[_i_can];
        _iPS = spac.plant_ps[_i_can];

        # set gsw to 0 or iterate for 30 times to find steady state solution
        if _df_dir + _df_dif < 10
            _iPS.APAR .= 0;
            _iPS.g_sw .= 0;
            gsw_control!(spac.photo_set, _iPS, _iEN);
        else
            for _ in 1:30
                gas_exchange!(spac.photo_set, _iPS, _iEN, GswDrive());
                update_gsw!(spac, spac.stomata_model, _i_can, FT(120); β = _βm);
                gsw_control!(spac.photo_set, _iPS, _iEN);
            end;
        end;
    end;

    # calculate the SIF if there is sunlight
    if _df_dir + _df_dif >= 10
        update_sif!(spac);
        dfr.SIF683 = SIF_WL(spac.can_rad, spac.wl_set, FT(682.5));
        dfr.SIF740 = SIF_740(spac.can_rad, spac.wl_set);
        dfr.SIF757 = SIF_WL(spac.can_rad, spac.wl_set, FT(758.7));
        dfr.SIF771 = SIF_WL(spac.can_rad, spac.wl_set, FT(770.0));
        dfr.NDVI   = NDVI(spac.can_rad, spac.wl_set);
        dfr.EVI    = EVI(spac.can_rad, spac.wl_set);
        dfr.NIRv   = NIRv(spac.can_rad, spac.wl_set);
    end;

    # save the total flux into the DataFrame
    dfr.F_H2O = T_VEG(spac);
    dfr.F_CO2 = CNPP(spac);
    dfr.F_GPP = GPP(spac);

    return nothing
end


"""

    run_model!(spac::SPACMono{FT}, df::DataFrame, nc_out::String) where {FT<:AbstractFloat}

Run CliMA Land at a site for the enture year, given
- `spac` Soil plant air continuum struct
- `df` Weather driver dataframe
- `nc_out` File path to save the model output

"""
function run_model!(spac::SPACMono{FT}, df::DataFrame, nc_out::String) where {FT<:AbstractFloat}
    _in_rad_bak = deepcopy(spac.in_rad);
    _in_dir     = _in_rad_bak.E_direct' * spac.wl_set.dWL / 1000;
    _in_dif     = _in_rad_bak.E_diffuse' * spac.wl_set.dWL / 1000;
    _deepcopies = [_in_rad_bak, _in_dir, _in_dif];
    _beta_g     = BetaGLinearPsoil{FT}();

    # set up memory
    _spac_mem = SPACMemory{FT}();

    # iterate through the time steps
    for _dfr in eachrow(df)
        prescribe_parameters!(spac, _dfr, _spac_mem, _deepcopies);
        run_time_step!(spac, _dfr, _beta_g);
    end;

    # save simulation results to hard drive
    save_nc!(nc_out, df[:, DF_VARIABLES]);

    return nothing
end


@time dict = load("debug.jld2");
@time wddf = prepare_wd(dict, "debug.nc");
@time spac = prepare_spac(dict);
@time run_model!(spac, wddf, "debug.output.nc");
