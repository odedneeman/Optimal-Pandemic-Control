function ret = getDLDx(x, beta, curTime)
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
    endoResponse = curKappa * delta * theta * R;
    maxEndo = 1 - minN;
    if endoResponse > maxEndo
        endoResponse = maxEndo;
    end

    dLDS = 0;
    dLDE = 0;
    dLDI = gumbelParam * exp(-r * curTime) * (-phi * (-log(n) - 0.2 * (1 - n^5/(1 - endoResponse)^5)...
            + log(1 - endoResponse)) + phi * (log(w) - 0.2));
    dLDR = gumbelParam * exp(-r * curTime) * chi * delta * theta + ...
        gumbelParam * exp(-r * curTime) * (1 - D - phi * I)*(n^5 * endoResponse/(1 - endoResponse)^6  ...
        -curKappa * delta * theta / (1 - endoResponse));
    dLDD = gumbelParam * exp(-r * curTime) * (log(n) + 0.2 * (1 - n^5/(1 - endoResponse)^5) ...
            -log(1 - endoResponse)+ log(w) - 0.2);
    dLDC = 0;
    ret = [dLDS, dLDE, dLDI, dLDR, dLDD, dLDC];

end