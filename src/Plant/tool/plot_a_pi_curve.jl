using LaTeXStrings
using PyPlot

# this function plot An-Pi curve
function plot_a_pi_curve(;
                         v25::Number = 80.0,
                         j25::Number = 135.0,
                      Γ_star::Number = 2.5,
                         tem::Number = 298.15,
                         par::Number = 1200.0,
                        p_O₂::Number = 21278.25,
                         r25::Number = false,
                        unit::Number = "K")
    list_p,list_an,list_ag = get_a_pi_curve(
                                            v25 = v25,
                                            j25 = j25,
                                         Γ_star = Γ_star,
                                            tem = tem,
                                            par = par,
                                           p_O₂ = p_O₂,
                                            r25 = r25,
                                           unit = unit)

    # plot the data
    clf()
    plot(list_p, list_ag, color=:red  , label=L"$A_\mathrm{gross}$")
    plot(list_p, list_an, color=:green, label=L"$A_\mathrm{net}$"  )
    xlabel(L"$P_\mathrm{i}$ (Pa)"                                   , fontsize=16)
    ylabel(L"$A_\mathrm{net}$ ($\mathrm{\mu}$mol m$^{-2}$ s$^{-1}$)", fontsize=16)
    legend(loc="lower right")
end
