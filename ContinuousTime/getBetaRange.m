function ret = getBetaRange(curTime, x)
    global maxN;
    global minN;
    global delta;
    global theta;
    S = x(1);
    E = x(2);
    I = x(3);
    R = x(4);
    D = x(5);
    C = x(6);
    curKappa = getTVKappa(curTime);
    curEndo = curKappa * delta * theta * R;
    maxNWithEndo = maxN - curEndo;
    if maxNWithEndo < minN
        %disp("Endo response error!");
        maxNWithEndo = minN;
    end
    maxBeta = getBetaFromN(maxNWithEndo, curTime);
    minBeta = getBetaFromN(minN, curTime);
    ret = [minBeta, maxBeta];
end