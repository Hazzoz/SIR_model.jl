using SIR_model
using Test

@testset "BasicSIR" begin
    @test solve_SIR(5000,1,0,60,BasicSIR(0.3,0.5)) == 1
end

@testset "SIRForceOfInfection" begin
    @test solve_SIR(5000,1,0,60,SIRForceOfInfection(0.3,0.5,10)) == 1
end

@testset "SIRHerdImmunity" begin
    @test solve_SIR(5000,1,0,60,SIRHerdImmunity(0.3,0.5,10)) == 1
end
