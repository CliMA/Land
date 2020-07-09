function Yujie111SetSoilType(node::Yujie111{FT}, soil_type) where {FT<:AbstractFloat}
    st = lowercase(soil_type)
    if st=="sand"
        p_ssat = 0.121*9.8*998.0*1E-6
        c_ssat = 0.395
        b_ssat = 4.05
        k_ssat = 0.634
        println("soil type Sand applied...")
    elseif (st=="sandy loam") | (st=="sandyloam")
        p_ssat = 0.218*9.8*998.0*1E-6
        c_ssat = 0.435
        b_ssat = 4.90
        k_ssat = 0.125
        println("soil type Sandy Loam applied...")
    elseif st=="loam"
        p_ssat = 0.478*9.8*998.0*1E-6
        c_ssat = 0.451
        b_ssat = 5.39
        k_ssat = 0.025
        println("soil type Loam applied...")
    elseif (st=="clay loam") | (st=="clayloam")
        p_ssat = 0.630*9.8*998.0*1E-6
        c_ssat = 0.476
        b_ssat = 8.52
        k_ssat = 0.009
        println("soil type Clay Loam applied...")
    elseif st=="clay"
        p_ssat = 0.405*9.8*998.0*1E-6
        c_ssat = 0.482
        b_ssat = 11.4
        k_ssat = 0.005
        println("soil type Clay applied...")
    else
        p_ssat = 0.630*9.8*998.0*1E-6
        c_ssat = 0.476
        b_ssat = 8.52
        k_ssat = 0.009
        println("soil type not supported, use default Clay Loam instead...")
    end
    node.p_ssat = p_ssat
    node.b_ssat = c_ssat
    node.c_ssat = b_ssat
    node.k_ssat = k_ssat
end
