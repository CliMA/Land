# include("./scripts/yujie_111_run_birch.jl")

# make birch node
node_birch = Yujie111Init()
node_birch.k_rhiz = 5.0E8
node_birch.k_root = 2228.0
node_birch.k_stem = 4926.7
node_birch.k_leaf = 5424.3
node_birch.b_root = 1.879
node_birch.c_root = 2.396
node_birch.b_stem = 2.238
node_birch.c_stem = 9.380
node_birch.b_leaf = 1.897
node_birch.c_leaf = 2.203
node_birch.k_sla  = 1.14
node_birch.h_soil = 0.0
node_birch.h_root = 0.0
node_birch.h_stem = 1.0
node_birch.h_leaf = 0.0
node_birch.laba   = 4758.5
node_birch.gaba   = 1000.0




# read dataset and run simulation
data_birch = CSV.read("./data/dataset_birch.csv", delim="\t")

function GetSimulation(k_rhiz=1E10, mae=true)
    obs_a = []
    obs_e = []
    obs_p = []
    mod_a = []
    mod_e = []
    mod_p = []
    row,col = size(data_birch)
    for i in 1:row
        # update the traits
        t_leaf = Float32(data_birch[i,2 ])
        vp_air = Float32(data_birch[i,3 ])
        p_co2  = Float32(data_birch[i,4 ])
        par    = Float32(data_birch[i,5 ])
        p_soil = Float32(data_birch[i,8 ])
        vmax   = Float32(data_birch[i,11])
        jmax   = Float32(data_birch[i,12])
        laba   = Float32(data_birch[i,10])
        k_tree = Float32(data_birch[i,17])
        tmp_node = deepcopy(node_birch)
        Yujie111UpdateLeaf(tmp_node, laba, vmax, jmax/vmax)
        Yujie111UpdateSoilFromP(tmp_node, p_soil)
        tmp_node.k_rhiz = k_rhiz
        tmp_node.k_root = k_tree * 1.863013458
        tmp_node.k_stem = k_tree * 4.11942096
        tmp_node.k_leaf = k_tree * 4.535504046

        # get optima
        envir = [par p_co2 vp_air t_leaf]
        e,p,a,c,g,t = Yujie111GetOptimalFKnowT(tmp_node, envir)

        # expand the comparison list
        if obs_a==[]
            obs_a = data_birch[i,14]
            obs_e = data_birch[i,18]
            obs_p = data_birch[i,16]
            mod_a = a
            mod_e = e
            mod_p = p
        else
            obs_a = [obs_a; data_birch[i,14]]
            obs_e = [obs_e; data_birch[i,18]]
            obs_p = [obs_p; data_birch[i,16]]
            mod_a = [mod_a; a]
            mod_e = [mod_e; e]
            mod_p = [mod_p; p]
        end
    end

    # calculate the mae
    obs_n = obs_p[obs_p.>0]
    mod_n = mod_p[obs_p.>0]
    std_a = mean(obs_a)
    std_e = mean(obs_e)
    std_p = mean(obs_n)
    mae_a = mean( abs.(obs_a-mod_a) / std_a )
    mae_e = mean( abs.(obs_e-mod_e) / std_e )
    mae_p = mean( abs.(obs_n-mod_n) / std_p )

    # return the results
    if mae
        return [mae_a mae_e mae_p]
    else
        return [obs_a obs_e obs_p mod_a mod_e mod_p]
    end
end

#=
list_rhiz = []
list_maes = []
for i in 1:100
    r_rhiz = i * 1E7
    ma,me,mp = GetSimulation(r_rhiz, true)
    maes = (ma + me + mp)/3.0
    println( @sprintf "%3d\t%8.3f\t%8.3f\t%8.3f\t%8.5f" i ma me mp maes )
    push!(list_rhiz, r_rhiz)
    push!(list_maes, maes  )
end
plot(list_rhiz, list_maes)
=#




# plot the figure
mat = GetSimulation(5.6E8, false)
p1 = scatter(mat[:,1], mat[:,4],color=:black,
             xlabel="obs A", ylabel="mod A",
             xlim=(0,28), ylim=(0,28));
p1 = plot!([0;28],[0;28],color=:black);
p2 = scatter(mat[:,2], mat[:,5],color=:black,
             xlabel="obs E", ylabel="mod E",
             xlim=(0,820), ylim=(0,820));
p2 = plot!([0;820],[0;820],color=:black);
obs_p = mat[:,3]
mod_p = mat[:,6]
mask  = (obs_p .> 0)
p3 = scatter(obs_p[mask], mod_p[mask],color=:black,
             xlabel="obs P", ylabel="mod P",
             xlim=(0.5,2.7), ylim=(0.5,2.7));
p3 = plot!([0.5;2.7],[0.5;2.7],color=:black);
fg = plot(p1,p2,p3,
          layout=(1,3), size=(1500,500), linecolor=:black, legend=false)
savefig(fg, "birch_best.pdf")

writedlm("./birch_simulation.txt", mat)
