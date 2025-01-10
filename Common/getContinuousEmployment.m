function ret = getContinuousEmployment(t1, n1, finalTime, simulationDt)
    finalStep = fix(finalTime / simulationDt);
    ret = zeros(1, finalStep + 1);
    curStage = 1;
    for curStep = 1 : finalStep + 1
        curTime = (curStep - 1) * simulationDt;
        if curStage <= length(t1) && curTime >= t1(curStage)
            curStage = curStage + 1;
        end
        ret(curStep) = n1(curStage);
    end
end

