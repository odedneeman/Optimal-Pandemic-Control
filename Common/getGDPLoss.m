function ret = getGDPLoss(xList, nList, simulationDt, finalTime)
    global phi;
    global r;
    global w;
    tList = 0:simulationDt:finalTime;
    expTerm = exp(-r * tList);
    stateTerm = (xList(5,:) + phi * xList(3,:)) + (1 - xList(5,:) - phi * xList(3,:)) .* (1 - nList);
    integrand = expTerm .* stateTerm;
    integral = sum(integrand * simulationDt);
    ret = (1/365) * integral;
end