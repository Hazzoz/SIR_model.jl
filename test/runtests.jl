using SIR_model
using Test

@testset "BasicSIR" begin
    #@test solve_SIR(5000,1,0,60,BasicSIR(0.3,0.5)) == 1
end

@testset "SIRForceOfInfection" begin
    #@test solve_SIR(5000,1,0,60,SIRForceOfInfection(0.3,0.5,10)) == 1
end

@testset "SIRHerdImmunity" begin
    sol = solve_SIR(5000,1,0,60,SIRHerdImmunity(0.3,0.5,10))

    R0 = 0.3*10/0.5
    pc = 1 - 1/R0

    check = 0
    limit = 0
    for i in sol.u
        if sum(i)*pc <= i[3] && check == 0
            limit = i[1]
            check = 1
        end
        if check != 0
            @test i[1] == limit
        end
    end
end
