function ret = costFunction(x, beta, curTime)
    global sigma;
    global gamma;
    global theta;
    global delta;
    global phi;
    global chi;
    global r;
    global w;
    global kappa;
    global minN;
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
    if n < minN
        n = minN;
    end
    endoResponse = curKappa * delta * theta * R;
    maxEndo = 1 - minN;
    if endoResponse > maxEndo
        endoResponse = maxEndo;
    end
    maxNWithEndo = 1 - endoResponse;
    if n > maxNWithEndo
        n = maxNWithEndo;
    end
    
    ret = gumbelParam * exp(-r * curTime) * ((1 - D - phi * I) * (-log(n) - 0.2 * (1 - n^5/(1 - endoResponse)^5)...
        + log(1 - endoResponse)) ...
        + (D + phi * I) * (log(w) - 0.2) + chi * delta * theta * R);

end