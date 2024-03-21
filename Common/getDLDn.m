function ret = getDLDn(x, beta, curTime, n)
    global sigma;
    global gamma;
    global theta;
    global delta;
    global phi;
    global chi;
    global r;
    global w;
    global kappa;
    curKappa = getTVKappa(curTime);
    gumbel = getGumbelCDF(curTime);
    gumbelParam = 1 - gumbel;
    S = x(1);
    E = x(2);
    I = x(3);
    R = x(4);
    D = x(5);
    C = x(6);
    %n = getNFromBeta(beta, curTime);
    endoResponse = curKappa * delta * theta * R;
    maxEndo = 1 - minN;
    if endoResponse > maxEndo
        endoResponse = maxEndo;
    end

    ret = gumbelParam * exp(-r * curTime) * ((1 - D - phi * I) * (-1/n - ...
        0.2 * ( - 5 * (n^4)/(1 - endoResponse)^5)));
end