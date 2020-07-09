# include("scripts/gnight/test_dade_drde.jl")

node = Yujie111Init()
e_crit = Yujie111GetECrit(node)

# to calculate drde
list_e = []
list_r = []
list_s = []
temp_e = 0.0
while true
    global e_crit,temp_e,list_e,list_r,list_s
    tlef    = Yujie111GetLeafTem(node, 25.0, 0.0, temp_e*0.0154321, false, 1.0)
    anet    = -GetRespiration(node.v_max, tlef)
    temp_e += 1.0
    if (temp_e>=e_crit) | (tlef<=0.0)
        break
    end
    tlef_de = Yujie111GetLeafTem(node, 25.0, 0.0, temp_e*0.0154321, false, 1.0)
    anet_de = -GetRespiration(node.v_max, tlef_de)
    drde    = anet_de - anet
    if list_e == []
        list_e = [temp_e]
        list_r = [anet  ]
        list_s = [drde  ]
    else
        list_e = [list_e; temp_e]
        list_r = [list_r; anet  ]
        list_s = [list_s; drde  ]
    end
end

# to calculate dade
envir  = [35.0, 40.0, 2.5, 3.0]
zenith = 30.0
r_all  = 1000.0
dade = Yujie111GetOptimaldAdE(node, envir, zenith, r_all)
println(dade)

comp_fact = 2.5
plot(list_e, list_s*comp_fact)
plot!(list_e, list_s*0.0 .+ dade)