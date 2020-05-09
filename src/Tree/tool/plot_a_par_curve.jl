using LaTeXStrings
using PyPlot

# this function plot An-Pi curve
function plot_a_par_curve(v25=80.0, j25=135.0, gamma=2.5, gsc=0.2, ca=40.0, tem=298.15, p_atm=101325.0,p_o2=21278.25, r25=false, unit="K")
    list_par,list_an,list_ag,list_pi = get_a_par_curve(v25,j25,gamma,gsc,ca,tem,p_atm,p_o2,r25,unit)

    # plot the data
    clf()
    plot(list_par, list_ag, color=:red  , label=L"$A_\mathrm{gross}$")
    plot(list_par, list_an, color=:green, label=L"$A_\mathrm{net}$"  )
    xlabel(L"PAR ($\mathrm{\mu}$mol m$^{-2}$ s$^{-1}$)"             , fontsize=16)
    ylabel(L"$A_\mathrm{net}$ ($\mathrm{\mu}$mol m$^{-2}$ s$^{-1}$)", fontsize=16)
    legend(loc="center right")

    twinx()
    plot(list_par, list_pi, color=:blue , label=L"$P_\mathrm{i}$"    )
    ylabel(L"$P_\mathrm{i}$ (Pa)", fontsize=16, color=:blue)
end
