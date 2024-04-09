function ret = hamiltonianDerivative(x, beta, p, curTime)
    global r;
    global phi;
    global alpha;
    global betaN;
    global kappa;
    global delta;
    global theta;
    curKappa = getTVKappa(curTime);
    gumbel = getGumbelCDF(curTime);
    gumbelParam = 1 - gumbel;
    S = x(1);
    E = x(2);
    I = x(3);
    R = x(4);
    D = x(5);
    C = x(6);

    n = getNFromBeta(beta, curTime);
    
    dLDN = gumbelParam * exp(-r * curTime) * ((1 - D - phi * I)*(-1/n + n^4/(1 - curKappa * delta * theta * R)^5));
    dBetaDN = betaN * alpha * ((1 - n)^(alpha - 1));

    dLDBeta = dLDN / dBetaDN;
    dDyanmicsDBeta = [-I * S; I * S; 0; 0; 0; 0];

    ret = dLDBeta + p' * dDyanmicsDBeta;

end