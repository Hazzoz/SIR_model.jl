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
#   - lambdas = recorded array of lambdas with timestamps
# - t = timespan
###############################################################
function SIR!(dP, P, params::BasicSIR, t)
    N = P[1] + P[2] + P[3] # Total Population
    lambda = P[2]*params.beta # Force of Infection
    dP[1] = -lambda/N*P[1] # Change in Susceptible Population
    dP[2] = lambda/N*P[1] - params.gamma*P[2] # Change in Infected Population
    dP[3] = params.gamma*P[2] # Change in Recovered Population

    # Store lambda at this time step
    push!(params.lambdas, (lambda, t)) # Storing lambda values
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
#   - lambdas = recorded array of lambdas with timestamps
# - t = timespan
###############################################################
function SIR!(dP, P, params::SIRForceOfInfection, t)
    N = P[1] + P[2] + P[3] # Total Population
    lambda = P[2]/N*params.beta*params.contacts # Force of Infection
    dP[1] = -lambda*P[1] # Change in Susceptible Population
    dP[2] = lambda*P[1] - params.gamma*P[2] # Change in Infected Population
    dP[3] = params.gamma*P[2] # Change in Recovered Population

    # Store lambda at this time step
    push!(params.lambdas, (lambda, t)) # Storing lambda values
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
#   - lambdas = recorded array of lambdas with timestamps
# - t = timespan
###############################################################
function SIR!(dP, P, params::SIRHerdImmunity, t)
    N = P[1] + P[2] + P[3] # Total Population
    R0 = params.contacts*params.beta/params.gamma
    pc = 1 - 1/R0 # Herd immunity threshold
    lambda = P[2]/N*params.beta*params.contacts  # Force of Infection
    #if P[3]+P[2] >= pc*N
    if P[3] >= pc*N # Check if recovered population exceeds herd immunity threshold
        lambda = 0
    end

    dP[1] = -lambda*P[1] # Change in Susceptible Population
    dP[2] = lambda*P[1] - params.gamma*P[2] # Change in Infected Population
    dP[3] = params.gamma*P[2] # Change in Recovered Population

    # Store lambda at this time step
    push!(params.lambdas, (lambda, t)) # Storing lambda values
end

###############################################################
# solve_SIR is a driver function that chooses the required SIR
# model
# Inputs:
# - S0 = Initial Susceptible Population
# - I0 = Initial Infected Population
# - R0 = Initial Recovered Population
# - days = No. of days modelled 
# - params = array of other necessary parameters
#   - beta = transmission chance of any interaction
#   - gamma = recovery rate
#   - contacts = number of daily contacts a person has
###############################################################
function solve_SIR(S0, I0, R0, days, params)
    P0 = [S0, I0, R0] # Initial populations vector

    tspan = (0, days) # Time span tuple

    solution = solve(ODEProblem(SIR!, P0, tspan, params)) # Solve the ODE with given parameters and timespan

    return solution, params.lambdas # Return values
end

###############################################################
# plot_SIR is a driver function to solve and plot the SIR model
# Inputs:
# - S0 = Initial Susceptible Population
# - I0 = Initial Infected Population
# - R0 = Initial Recovered Population
# - days = No. of days modelled 
# - params = array of other necessary parameters
#   - beta = transmission chance of any interaction
#   - gamma = recovery rate
#   - contacts = number of daily contacts a person has
###############################################################
function plot_SIR(S0, I0, R0, days, params)
    solution, lambdas = solve_SIR(S0, I0, R0, days, params) # Solve the SIR model

    plot(solution, xlabel="Time", ylabel="Population", title="Solution", labels=["Susceptible" "Infected" "Recovered"]) # Plot the model
end