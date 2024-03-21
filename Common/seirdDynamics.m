function ret = seirdDynamics(x, beta)
    global sigma;
    global gamma;
    global theta;
    global delta;
    % S - Susceptible; E - Exposed; I - Infected; R - Resolving; D - Death; 
    % C - Recovered.
    S = x(1);
    E = x(2);
    I = x(3);
    R = x(4);
    D = x(5);
    C = x(6);
    dS = - beta * I * S;
    dE = beta * I * S - sigma * E;
    dI = sigma * E - gamma * I;
    dR = gamma * I - theta * R;
    dD = delta * theta * R;
    dC = (1 - delta) * theta * R;
    
    ret = [dS; dE; dI; dR; dD; dC];
end

