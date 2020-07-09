# include("./scripts/leaf investment/yujie_111_optimal_leaf_investment.jl")
# make optimal node
ini_laba = 2000.0
ini_vmax = 60.0

opt_node = Yujie111Init()
Yujie111UpdateSoilFromSWC(opt_node, 1.0)

weat= DataFrame!(CSV.File("./data/2070.csv"))
opt_laba,opt_vmax,opt_prof = Yujie111GetOptimalInvestment(opt_node, weat)

Yujie111UpdateLeaf(opt_node, opt_laba, opt_vmax)
