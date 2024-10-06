using SIR_model
using Test

@testset "BasicSIR" begin
    # Solve the SIR model
    sol = solve_SIR(5000, 1, 0, 60, BasicSIR(0.3, 0.1))

    # Get the initial total population
    initial_total = sum(sol[1])

    # Check that the total population remains the same at every time step
    for i in sol
        @test sum(i) ≈ initial_total
    end
end

@testset "SIRForceOfInfection" begin
    #@test solve_SIR(5000,1,0,60,SIRForceOfInfection(0.3,0.5,10)) == 1
end

@testset "SIRHerdImmunity" begin
    #@test solve_SIR(5000,1,0,60,SIRHerdImmunity(0.3,0.5,10,0.8)) == 1
end
