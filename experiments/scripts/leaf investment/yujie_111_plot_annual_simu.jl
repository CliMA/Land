# include("experiments/scripts/leaf investment/yujie_111_plot_annual_simu.jl")

using BenchmarkTools
using CSV
using DataFrames
using Plants
using Photosynthesis

FT = Float64

# make node for year 2005
node = Plants.SPACSimple{FT}();
node.d_lati =  35.198284;
node.d_long = -111.651299;
node.d_alti =  2105.0;

c3_set = C3CLM(FT);

#=
factor = 0.2
node.b_root *= factor
node.b_stem *= factor
node.b_leaf *= factor
Yujie111UpdateLeaf(node, 1318.0, 77.9)
=#

#=
factor = 0.6
node.b_root *= factor
node.b_stem *= factor
node.b_leaf *= factor
Yujie111UpdateLeaf(node, 1612.0, 77.7)
=#

#=
factor = 1.5
node.b_root *= factor
node.b_stem *= factor
node.b_leaf *= factor
#Yujie111UpdateLeaf(node, 1701.0, 77.7)
=#

#=
factor = 2.0
node.b_root *= factor
node.b_stem *= factor
node.b_leaf *= factor
#Yujie111UpdateLeaf(node, 1733.0, 76.4)
=#

#= 0.2X
node.k_root *= 0.2
node.k_stem *= 0.2
node.k_sla  *= 0.2
Yujie111UpdateLeaf(node, 1509.0, 74.6)
=#

#=
node.k_root *= 0.6
node.k_stem *= 0.6
node.k_sla  *= 0.6
Yujie111UpdateLeaf(node, 1509.0, 81.2)
=#

#=
node.k_root *= 1.5
node.k_stem *= 1.5
node.k_sla  *= 1.5
Yujie111UpdateLeaf(node, 1220.0, 78.8)
# =#

#= 2.0X
node.k_root *= 2.0
node.k_stem *= 2.0
node.k_sla  *= 2.0
Yujie111UpdateLeaf(node, 1627.0, 79.3)
=#

# #= default
Plants.Yujie111UpdateLeaf(node, c3_set, 1193.0, 88.9);
# =#

Plants.soil_moisture_swc!(node, 1.0);

#weat_hist = DataFrame!(CSV.File("./data/Flagstaff_histo_sim0_GS.csv"))
weat_hist = DataFrame!(CSV.File("/home/wyujie/Data/USAForest20Sites/Flagstaff/Flagstaff_histo_sim0_GS.csv"))
weat_mask = (weat_hist.Year .== 2005);
weat_2005 = weat_hist[weat_mask,:];

tmp_node = deepcopy(node);

# initialize data frame to store the results
df = DataFrame()
df[!, "Time"  ]  = weat_2005.Day + weat_2005.Hour / FT(24)
df[!, "T_air" ] .= FT(0)
df[!, "D_air" ] .= FT(0)
df[!, "Wind"  ] .= FT(0)
df[!, "Rain"  ] .= FT(0)
df[!, "Ca"    ] .= FT(0)
df[!, "SWC"   ] .= FT(0)
df[!, "P_soil"] .= FT(0)
df[!, "H_sun" ] .= FT(0)
df[!, "A_net" ] .= FT(0)
df[!, "LAI_sl"] .= FT(0)
df[!, "PAR_sl"] .= FT(0)
df[!, "RAD_sl"] .= FT(0)
df[!, "E_sl"  ] .= FT(0)
df[!, "P_sl"  ] .= FT(0)
df[!, "An_sl" ] .= FT(0)
df[!, "Ag_sl" ] .= FT(0)
df[!, "C_sl"  ] .= FT(0)
df[!, "G_sl"  ] .= FT(0)
df[!, "T_sl"  ] .= FT(0)
df[!, "LAI_sh"] .= FT(0)
df[!, "PAR_sh"] .= FT(0)
df[!, "RAD_sh"] .= FT(0)
df[!, "E_sh"  ] .= FT(0)
df[!, "P_sh"  ] .= FT(0)
df[!, "An_sh" ] .= FT(0)
df[!, "Ag_sh" ] .= FT(0)
df[!, "C_sh"  ] .= FT(0)
df[!, "G_sh"  ] .= FT(0)
df[!, "T_sh"  ] .= FT(0)


@time Plants.simulate_growing_season!(tmp_node, c3_set, weat_2005, df);

tmp_node = deepcopy(node);
@time Plants.simulate_growing_season!(tmp_node, c3_set, weat_2005, df);

# save matrix
# writedlm("./annual_simulation.txt", mat)

#=
1. time
2. T_air
3. D_air
4. Wind
5. Rain
6. soil water content
7. P_soil
8. solar height angle
9. A_net

10. LAI_sublit
11. PAR_sunlit
12. E_sunlit
13. P_sunlit
14. A_sunlit
15. C_sunlit
16. G_sunlit
17. T_sunlit

18. LAI_shade
19. PAR_shade
20. E_shade
21. P_shade
22. A_shade
23. C_shade
24. G_shade
25. T_shade
=#

#=
p2  = plot(mat[:,1], mat[:,2], ylabel=L"$T_\mathrm{air}$" , xlim=(130,280))
p3  = plot(mat[:,1], mat[:,3], ylabel=L"$D_\mathrm{air}$" , xlim=(130,280))
p4  = plot(mat[:,1], mat[:,4], ylabel="Wind"              , xlim=(130,280))
p5  = plot(mat[:,1], mat[:,5], ylabel="Precipitation"     , xlim=(130,280))
p6  = plot(mat[:,1], mat[:,6], ylabel="SWC"               , xlim=(130,280))
p7  = plot(mat[:,1], mat[:,7], ylabel=L"$P_\mathrm{soil}$", xlim=(130,280))
p8  = plot(mat[:,1], mat[:,8], ylabel="Solar height"      , xlim=(130,280))
p9  = plot(mat[:,1], mat[:,9], ylabel=L"A_\mathrm{net}"   ,
                               xlabel="Day in Year"       , xlim=(130,280))
fg1 = plot(p2,p3,p4,p5,p6,p7,p8,p9,
           layout=(8,1), size=(1500,1200), linecolor=:black, legend=false)
savefig(fg1, "season_1.pdf")

p10 = plot(mat[:,1], mat[:,10], ylabel=L"\mathrm{LAI_{sunlit}}", xlim=(130,280))
p11 = plot(mat[:,1], mat[:,11], ylabel=L"\mathrm{PAR_{sunlit}}", xlim=(130,280))
p12 = plot(mat[:,1], mat[:,12], ylabel=L"E_\mathrm{sunlit}"    , xlim=(130,280))
p13 = plot(mat[:,1], mat[:,13], ylabel=L"P_\mathrm{sunlit}"    , xlim=(130,280))
p14 = plot(mat[:,1], mat[:,14], ylabel=L"A_\mathrm{sunlit}"    , xlim=(130,280))
p15 = plot(mat[:,1], mat[:,15], ylabel=L"C_\mathrm{sunlit}"    , xlim=(130,280))
p16 = plot(mat[:,1], mat[:,16], ylabel=L"G_\mathrm{sunlit}"    , xlim=(130,280))
p17 = plot(mat[:,1], mat[:,17], ylabel=L"T_\mathrm{sunlit}"    ,
                                xlabel="Day in Year"           , xlim=(130,280))
fg2 = plot(p10,p11,p12,p13,p14,p15,p16,p17,
           layout=(8,1), size=(1500,1200), linecolor=:black, legend=false)
savefig(fg2, "season_2.pdf")

p18 = plot(mat[:,1], mat[:,18], ylabel=L"\mathrm{LAI_{shade}}", xlim=(130,280))
p19 = plot(mat[:,1], mat[:,19], ylabel=L"\mathrm{PAR_{shade}}", xlim=(130,280))
p20 = plot(mat[:,1], mat[:,20], ylabel=L"E_\mathrm{shade}"    , xlim=(130,280))
p21 = plot(mat[:,1], mat[:,21], ylabel=L"P_\mathrm{shade}"    , xlim=(130,280))
p22 = plot(mat[:,1], mat[:,22], ylabel=L"A_\mathrm{shade}"    , xlim=(130,280))
p23 = plot(mat[:,1], mat[:,23], ylabel=L"C_\mathrm{shade}"    , xlim=(130,280))
p24 = plot(mat[:,1], mat[:,24], ylabel=L"G_\mathrm{shade}"    , xlim=(130,280))
p25 = plot(mat[:,1], mat[:,25], ylabel=L"T_\mathrm{shade}"    ,
                                xlabel="Day in Year"          , xlim=(130,280))
fg3 = plot(p18,p19,p20,p21,p22,p23,p24,p25,
           layout=(8,1), size=(1500,1200), linecolor=:black, legend=false)
savefig(fg3, "season_3.pdf")
=#
