function ret = getContinuousInputFromDiscrete(n1, t1, simulationDt, finalStep)
    continuousNList = zeros(1, finalStep);
    curN = 1;
    curT = 1;
    for curStep = 1:finalStep
        curTime = (curStep - 1)* simulationDt;
        if curT <= length(t1)
            if curTime >= t1(curT)
                curT = curT + 1;
                curN = curN + 1;
            end
        end
        continuousNList(curStep) = n1(curN);
    end
    ret = continuousNList;
end

