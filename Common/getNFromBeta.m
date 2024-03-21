function ret = getNFromBeta(beta, curTime)
    % Assuming n_ss = 1
    global betaW;
    global betaN;
    global betaLambda;
    global lambda;
    global alpha;
    ret = 1 - ((beta - betaW - betaLambda * exp( -lambda * curTime ))/(-betaN))^(1/alpha); 
end