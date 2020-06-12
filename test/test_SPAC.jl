# Diurnal cycle
# TODO more functions to be added
@testset "CanopyRT --- Diurnal cycle" begin
    for FT in [Float32, Float64]
        wl_set    = RT.create_wl_para_set(FT)
        canopy_rt = RT.Canopy4RT{FT, 20, 3.0}()
        canRad_rt = RT.CanopyRadiation{FT, wl_set.nwl, wl_set.nWlF, length(canopy_rt.litab), length(canopy_rt.lazitab), canopy_rt.nlayers}()
        canOpt_rt = RT.create_canopy_optical(FT, wl_set.nwl, canopy_rt.nlayers, length(canopy_rt.lazitab), length(canopy_rt.litab); using_marray=false)
        sunRad_rt = RT.create_incoming_radiation(FT, wl_set.swl)
        soil      = RT.SoilOpti{FT}(wl_set.wl, FT(0.2)*ones(FT, length(wl_set.wl)), FT[0.1], FT(290.0))
        angles    = RT.SolarAngles{FT}()

        arrayOfLeaves = [RT.create_leaf_bio(FT, wl_set.nwl, wl_set.nWlE, wl_set.nWlF) for i in 1:canopy_rt.nlayers]
        for i in 1:canopy_rt.nlayers
            RT.fluspect!(arrayOfLeaves[i], canopy_rt, wl_set)
        end

        RT.compute_canopy_geometry!(canopy_rt, angles, canOpt_rt)
        RT.compute_canopy_matrices!(arrayOfLeaves, canOpt_rt)
        RT.simulate_short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil)
        RT.derive_canopy_fluxes!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil, arrayOfLeaves, wl_set)

        tree = PT.Tree{FT,5,20,325}()
        PT.update_tree_e_crit!(tree)

        SPAC.update_canopy_from_rt_module!(tree, canopy_rt, canOpt_rt, canRad_rt)
        canopy_par_0 = []
        for layer in 1:tree.n_canopy
            push!(canopy_par_0, tree.canopy_list[layer].par_list[:])
        end

        ## A diurnal cycle for radiation and Tair
        Deltat = FT(60)
        Tmean  = FT(295.15)
        DeltaT = 3
        omega  = 2 * π / (24*3600)
        t      = range(0, stop=24*3600, step=Deltat)
        phi_t  = omega * t - π * ones(FT, size(t)) / 2
        PARr_t = zeros(FT, size(t))
        Tair_t = zeros(FT, size(t))
        N      = length(PARr_t)

        for i = 1:N
            PARr_t[i] = max( sin(phi_t[i]), FT(0.0) )
            Tair_t[i] = Tmean + DeltaT * sin(phi_t[i] - π/3)
        end

        # TODO change zenith angle and PAR accordingly
        # TODO add leaf energy balance
        # TODO so slow here? switch to p_i orriented later
        for i in [collect(1:20); collect(500:520)]
            for layer in 1:tree.n_canopy
                tree.canopy_list[layer].par_list .= canopy_par_0[layer] .* PARr_t[i]
            end
            PT.update_tree_with_time!(tree, Deltat, tree.stomata_scheme)
            #println( i, "\t", PARr_t[i], "\t", mean(tree.canopy_list[end].gsw_list) )
        end

        @test true
    end
end
