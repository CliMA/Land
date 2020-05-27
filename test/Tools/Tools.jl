using LaTeXStrings
using Parameters
using Printf
using PyPlot

using Land.Plant

@unpack FT,
        Tree,
        get_a_par_curve,
        get_a_pi_curve = Plant

include("plot_a_par_curve.jl"     )
include("plot_a_pi_curve.jl"      )
include("visualize_struct_tree.jl")
