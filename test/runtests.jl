using SIR_model
using Test

@testset "BasicSIR" begin
    # Solve the SIR model
    sol, lambdas = solve_SIR(5000, 1, 0, 60, BasicSIR(0.3, 0.1, []))

    # Get the initial total population
    initial_total = sum(sol.u[1])

    # Check that the total population remains the same at every time step
    for i in sol.u
        @test sum(i) ≈ initial_total
    end
end

@testset "SIRForceOfInfection" begin
    # Solve the SIR model and obtain the lambda values
    sol, lambdas = solve_SIR(5000,1,0,60,SIRForceOfInfection(0.3,0.5,10,[]))

    # Loop through the time steps and check λ calculation
    for i in 1:length(sol.t)
        t = sol.t[i]
        lambda = 0
        for j in lambdas
            if j[2] == t
                lambda = j[1]
                break
            end
        end

        # Manually compute λ as β * I / N
        N = 5001  # Total population (constant)
        lambda_calc = 0.3 * 10 * sol.u[i][2] / N

        # Test if the manually computed λ matches the model’s λ
        @test lambda ≈ lambda_calc atol = 0.2
    end
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
