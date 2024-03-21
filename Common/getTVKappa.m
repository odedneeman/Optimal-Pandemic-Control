function ret = getTVKappa(curTime)
    global kappa;
    global phiKappa;
    global muF;
    global sigmaF;
    ret = kappa * (1 - (1 - phiKappa) * normcdf((curTime - muF)/sigmaF));
end