using SIR_model
using Test

@testset "BasicSIR" begin
    @test solve_SIR(5000,1,0,10,0.3,0.5,60) == 1
end

@testset "SIRForceOfInfection" begin
    @test solve_SIR(5000,1,0,10,0.3,0.5,60) == 1
end

@testset "SIRHerdImmunity" begin
    @test solve_SIR(5000,1,0,10,0.3,0.5,60) == 1
end
