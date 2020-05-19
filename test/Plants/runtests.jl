
using Land.Plant

@testset "Plants" begin
      
      FT = Float32
      v25::FT = FT(80.0)
      j25::FT = FT(135.0)
   Γ_star::FT = FT(2.5)
      gsc::FT = FT(0.1)
      p_a::FT = FT(40.0)
      tem::FT = FT(298.15)
    p_atm::FT = FT(101325.0)
     p_O₂::FT = FT(21278.25)
      r25::FT = FT(Inf)
    # breaks:
    # list_par,list_an,list_ag,list_pi = Plant.get_a_par_curve(
    #                                                          v25 = v25,
    #                                                          j25 = j25,
    #                                                       Γ_star = Γ_star,
    #                                                          gsc = gsc,
    #                                                          p_a = p_a,
    #                                                          tem = tem,
    #                                                        p_atm = p_atm,
    #                                                         p_O₂ = p_O₂,
    #                                                          r25 = r25)

end
