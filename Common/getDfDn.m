function ret = getDfDn(x, beta, curTime, n)
    global sigma;
    global gamma;
    global theta;
    global delta;
    
    S = x(1);
    E = x(2);
    I = x(3);
    R = x(4);
    D = x(5);
    C = x(6);
    %n = getNFromBeta(beta, curTime);
    ret = [-S * I;
            S * I;
            0;
            0;
            0;
            0];
    global alpha;
    global betaN;
    dBetaDn = alpha * betaN * (1 - n)^(alpha - 1);
    ret = ret .* dBetaDn;
end