module SIR_model

# Necessary Packages
Pkg.add("DifferentialEquations")
using DifferentialEquations

Pkg.add("Plots")
using Plots

# SIR Functions
include("mycode.jl")
export solve_SIR

end
