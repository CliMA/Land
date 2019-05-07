push!(LOAD_PATH, "../src/Utils/");  push!(LOAD_PATH, "../src/Leaf/")

using MathToolsMod, Leaf, LeafPhotosynthesisMod, Plots

tmax = 10000 # Time in seconds here
l = leaf_params{Float32}()
f = LeafPhotosynthesisMod.fluxes{Float32}()
l.dynamic_state = true
l.gstyp = 1

uu = zeros(tmax,10)
uu[1,2] = 0.07
uu[1,1] = 0.5
f.APAR = 100
l.rdleaf = 0.0
for c = 2:1:tmax
   #println(c)a
    if c==1000 f.APAR = 1500 end
    if c==5000 f.APAR = 100 end
       l.gs = uu[c-1,2]
       l.Kn = uu[c-1,1]

       uu[c-1,3]=l.ϕs
       uu[c-1,4]=l.Ci
       uu[c-1,5]=f.φ
       uu[c-1,6]=l.Kp
       uu[c-1,8]=f.an
       uu[c-1,7]=f.APAR
       #l.Kp = 4
       LeafPhotosynthesisMod.LeafPhotosynthesis(f,l,298.0)
       uu[c,1] = uu[c-1,1]+(l.Kn_ss-l.Kn)/10*1/60
       uu[c,2] = uu[c-1,2]+(l.gs_ss-l.gs)/15*1/60
       #println(l.Kn_ss, " ",  l.Ci, " ", f.φ)
end
