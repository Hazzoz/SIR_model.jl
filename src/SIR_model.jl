module SIR_model

    # Necessary Packages
    using DifferentialEquations
    using Plots

    # SIR Functions
    include("mycode.jl")
    export solve_SIR
    export BasicSIR, SIRForceOfInfection, SIRHerdImmunity

end
