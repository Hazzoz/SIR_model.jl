struct BasicSIR
    beta::Float64  # Transmission rate
    gamma::Float64  # Recovery rate
end

struct SIRForceOfInfection
    beta::Float64  # Transmission rate
    gammma::Float64  # Recovery rate
    contacts::Float64  # Force of infection
end

struct SIRHerdImmunity
    beta::Float64  # Transmission rate
    gammma::Float64  # Recovery rate
    contacts::Float64  # Force of infection
    herd::Float64  # Herd immunity threshold
end

###############################################################
# SIR! a function that represents to differential equations
# that define the SIR model
# Inputs:
# - dP = array of gradients for each population
# - P = array of current values for each population
# - params = array of other necessary parameters
#   - beta = transmission rate
#   - gamma = recovery rate
#   - 
# - t = timespan
###############################################################
function SIR!(dP, P, params::BasicSIR, t)
    beta, gamma = params
    N = P[1] + P[2] + P[3]
    lambda = P[2]*beta
    dP[1] = -lambda/N*P[1]
    dP[2] = lambda/N*P[1] - gamma*P[2]
    dP[3] = gamma*P[2]
end


###############################################################
# SIR! a function that represents to differential equations
# that define the SIR model
# Inputs:
# - dP = array of gradients for each population
# - P = array of current values for each population
# - params = array of other necessary parameters
# - t = timespan
###############################################################
function SIR!(dP, P, params::SIRForceOfInfection, t)
    beta, gamma, c = params
    N = P[1] + P[2] + P[3]
    lambda = P[2]/N*beta*c
    dP[1] = -lambda*P[1]
    dP[2] = lambda*P[1] - gamma*P[2]
    dP[3] = gamma*P[2]
end

###############################################################
# SIR! a function that represents to differential equations
# that define the SIR model
# Inputs:
# - dP = array of gradients for each population
# - P = array of current values for each population
# - params = array of other necessary parameters
# - t = timespan
###############################################################
function SIR!(dP, P, params::SIRHerdImmunity, t)
    c, beta, gamma, H = params
    N = P[1] + P[2] + P[3]
    lambda = P[2]/N*beta*c
    dP[1] = -lambda*P[1]
    dP[2] = lambda*P[1] - gamma*P[2]
    dP[3] = gamma*P[2]
end

###############################################################
# solve_SIR is a function that runs the ODE solver on the SIR!
# function using the parameters input into this command
# Inputs:
# - +Int: S0 = Initial Susceptible Population
# - +Int: I0 = Initial Infected Population
# - +Int: R0 = Initial Recovered Population
# - +Int: contacts = No. of daily contacts any person will have
# - +Float between 0-1: beta = Contraction chance of disease
# - +Float between 0-1: gamma = 
# - +Int: days = No. of days modelled 
###############################################################
function solve_SIR(S0, I0, R0, beta, gamma, days, contacts, herd)
    # Initial populations vector
    P0 = [S0, I0, R0]

    # Time span tuple
    tspan = (0, days)

    # Params
    if contacts 
        params::SIRHerdImmunity
        params.contacts = contacts
        params.herd = herd
    elseif contacts 
        params::SIRForceOfInfection
        params.contacts = contacts
    else
        params::BasicSIR
    end

    params.beta = beta
    params.gamma = gamma

    # Solve the ODE with given parameters and timespan
    solution = solve(ODEProblem(SIR!,P0,tspan,params))

    println(solution.u)

    # Plot the model
    #plot(solution, xlabel="Time", ylabel="Population", title="Solution", labels=["Susceptible" "Infected" "Recovered"])

end