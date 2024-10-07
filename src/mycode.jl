# Structs for different versions of the SIR models
struct BasicSIR
    beta::Float64  # Transmission rate
    gamma::Float64  # Recovery rate
    lambdas::Vector{Tuple{Float64,Float64}} # Recorded array of lambdas
end

struct SIRForceOfInfection
    beta::Float64  # Transmission chance of interaction
    gamma::Float64  # Recovery rate
    contacts::Float64  # Number of Daily Contacts
    lambdas::Vector{Tuple{Float64,Float64}} # Recorded array of lambdas
end

struct SIRHerdImmunity
    beta::Float64  # Transmission chance of interaction
    gamma::Float64  # Recovery rate
    contacts::Float64  # Number of Daily Contacts
    herd::Float64  # Herd immunity threshold
    lambdas::Vector{Tuple{Float64,Float64}} # Recorded array of lambdas
end

###############################################################
# This SIR! is a function that represents the differential
# equations that defines the a SIR model
# Inputs:
# - dP = array of gradients for each population
# - P = array of current values for each population
# - params = array of other necessary parameters
#   - beta = transmission rate
#   - gamma = recovery rate
# - t = timespan
###############################################################
function SIR!(dP, P, params::BasicSIR, t)
    N = P[1] + P[2] + P[3]
    lambda = P[2]*params.beta
    dP[1] = -lambda/N*P[1]
    dP[2] = lambda/N*P[1] - params.gamma*P[2]
    dP[3] = params.gamma*P[2]

    # Store lambda at this time step
    push!(params.lambdas, (lambda, t))
end

###############################################################
# This SIR! is a function that represents the differential
# equations that defines a more detailed force of infection
# version of the SIR model
# Inputs:
# - dP = array of gradients for each population
# - P = array of current values for each population
# - params = array of other necessary parameters
#   - beta = transmission chance of any interaction
#   - gamma = recovery rate
#   - contacts = number of daily contacts a person has
# - t = timespan
###############################################################
function SIR!(dP, P, params::SIRForceOfInfection, t)
    N = P[1] + P[2] + P[3]
    lambda = P[2]/N*params.beta*params.contacts
    dP[1] = -lambda*P[1]
    dP[2] = lambda*P[1] - params.gamma*P[2]
    dP[3] = params.gamma*P[2]

    # Store lambda at this time step
    push!(params.lambdas, (lambda, t))
end

###############################################################
# This SIR! is a function that represents the differential
# equations that defines an SIR model with a herd immunity
# threshold
# Inputs:
# - dP = array of gradients for each population
# - P = array of current values for each population
# - params = array of other necessary parameters
#   - beta = transmission chance of any interaction
#   - gamma = recovery rate
#   - contacts = number of daily contacts a person has
#   - herd = herd immunity threshold
# - t = timespan
###############################################################
function SIR!(dP, P, params::SIRHerdImmunity, t)
    N = P[1] + P[2] + P[3]
    R0 = params.contacts*params.beta/params.gamma
    pc = 1 - 1/R0
    lambda = P[2]/N*params.beta*params.contacts
    if P[3] >= pc*N
        lambda = 0
    end
    dP[1] = -lambda*P[1]
    dP[2] = lambda*P[1] - params.gamma*P[2]
    dP[3] = params.gamma*P[2]

    # Store lambda at this time step
    push!(params.lambdas, (lambda, t))
end

###############################################################
# solve_SIR is a driver function that chooses the required SIR
# model
# Inputs:
# - +Int: S0 = Initial Susceptible Population
# - +Int: I0 = Initial Infected Population
# - +Int: R0 = Initial Recovered Population
# - +Int: days = No. of days modelled 
# - params = array of other necessary parameters
#   - beta = transmission chance of any interaction
#   - gamma = recovery rate
#   - contacts = number of daily contacts a person has
#   - herd = herd immunity threshold
###############################################################
function solve_SIR(S0, I0, R0, days, params)
    # Initial populations vector
    P0 = [S0, I0, R0]

    # Time span tuple
    tspan = (0, days)

    # Solve the ODE with given parameters and timespan
    solution = solve(ODEProblem(SIR!, P0, tspan, params))

    return solution, params.lambdas
end

###############################################################
# plot_SIR is a driver function to solve and plot the SIR model
# Inputs:
# - +Int: S0 = Initial Susceptible Population
# - +Int: I0 = Initial Infected Population
# - +Int: R0 = Initial Recovered Population
# - +Float between 0-1: beta = Transmission rate/chance of interaction
# - +Float between 0-1: gamma = Recovery Rate
# - +Int: days = No. of days modelled 
# - +Int: contacts = No. of daily contacts any person will have
# - +Float between 0-1: herd = Herd immunity threshold
###############################################################
function plot_SIR(S0, I0, R0, days, params)
    # Solve the SIR model
    solution = solve_SIR(S0, I0, R0, days, params)

    # Plot the model
    plot(solution, xlabel="Time", ylabel="Population", title="Solution", labels=["Susceptible" "Infected" "Recovered"])
end