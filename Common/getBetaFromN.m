function ret = getBetaFromN(n, curTime)
    % Assuming n_ss = 1
    global betaW;
    global betaN;
    global alpha;
    global betaLambda;
    global lambda;
    ret = betaW - betaN * (1 - n)^alpha + betaLambda * exp(-lambda * curTime);
end