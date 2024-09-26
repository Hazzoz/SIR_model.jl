###############################################################
# SIR! a function that represents to differential equations
# that define the SIR model
# Inputs:
# - dP = array of gradients for each population
# - P = array of current values for each population
# - params = array of other necessary parameters
# - t = timespan
###############################################################
function SIR!(dP, P, params, t)
    c, beta, gamma, N = params
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
function solve_SIR(S0, I0, R0, contacts, beta, gamma, days)
    # Initial populations vector
    P0 = [S0, I0, R0]

    # Parameters vector (P[4] = total Population)
    params = [contacts, beta, gamma, S0 + I0 + R0]

    # Time span tuple
    tspan = (0, days)

    # Solve the ODE with given parameters and timespan
    solution = solve(ODEProblem(SIR,P0,tspan,params))

    # Plot the model
    plot(solution, xlabel="Time", ylabel="Population", title="Solution", labels=["Susceptible" "Infected" "Recovered"])
end