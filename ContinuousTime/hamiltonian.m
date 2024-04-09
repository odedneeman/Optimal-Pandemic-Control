function ret = hamiltonian(x, beta, p, curTime)
    costValue = costFunction(x, beta, curTime);
    costateValue = p' * seirdDynamics(x, beta);
    ret =  costValue + costateValue;
end