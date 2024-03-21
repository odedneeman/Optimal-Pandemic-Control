function ret = hamiltonianDoubleDerivative(x, beta, p, curTime)
    global r;
    global phi;
    global alpha;
    global betaN;
    global kappa;
    global delta;
    global theta;
    gumbel = getGumbelCDF(curTime);
    gumbelParam = 1 - gumbel;
    curKappa = getTVKappa(curTime);
    S = x(1);
    E = x(2);
    I = x(3);
    R = x(4);
    D = x(5);
    C = x(6);
    n = getNFromBeta(beta, curTime);
    dBetaDN = betaN * alpha * ((1 - n)^(alpha - 1));
    doubleDLDBeta = gumbelParam * (exp(-r * curTime) * (1 - D - phi * I) * (1/(n^2) + 4 * (n^3)/(1 - curKappa * delta * theta * R)^5) * dBetaDN + ...
        exp(-r * curTime) * ((1 - D - phi * I)*(-1/n + n^4/(1 - curKappa * delta * theta * R)^5)) * betaN * alpha * (alpha - 1) * ((1 - n)^(alpha - 2))) / ...
        ((dBetaDN)^3);
    ret = doubleDLDBeta;
end