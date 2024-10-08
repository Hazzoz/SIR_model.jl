module SIR_model

    # Necessary Packages
    using DifferentialEquations
    using Plots

    # SIR Functions
    include("mycode.jl")
    export solve_SIR, plot_SIR, get_vals, compare_vals
    export BasicSIR, SIRForceOfInfection, SIRHerdImmunity

end
