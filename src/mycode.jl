function SIR!(dP, P, params, t)
    c, beta, gamma, N = params
    lambda = P[2]/N*beta*c
    dP[1] = -lambda*P[1]
    dP[2] = lambda*P[1] - gamma*P[2]
    dP[3] = gamma*P[2]
end