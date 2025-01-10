function ret = getEquivalentLossAllSteps(xList, betaList, simulationDt, finalStep)
    global r;
    % finalStep is different in descrete system, for this function it needs
    % to be adjusted.
    finalStepAdj = finalStep;
    if mod(finalStep, 10) == 0
        finalStepAdj = finalStep + 1;
    end
    wValueList = zeros(1, finalStepAdj);
    for curStep = 1 : finalStepAdj
        curTime = (curStep - 1) * simulationDt;
        wValueList(curStep) = costFunctionForLoss(xList(:,curStep), ...
                            betaList(curStep), curTime);
    end
    equivNList = zeros(1, finalStepAdj);
    options = optimset('Display','off');
    for curStep = 1 : finalStepAdj
        nEquation = @(n)(-log(n) - 1/5*(1 - n.^5) - wValueList(curStep));
        curNValue = fsolve(nEquation, 0.5, options);
        equivNList(curStep) = curNValue;
    end

    % plot(0:0.01:635, equivNList)
    % hold on
    % plot(0:0.01:635,nList)
    % legend("equivalent employment rate")
    % legend(["equivalent employment rate", "true employment rate"])
    finalTime = (finalStepAdj - 1)*simulationDt;
    accumLoss = zeros(1, finalStepAdj);
    totalLoss = 0;
    for curStep = 1 : finalStepAdj
        curTime = (curStep - 1)*simulationDt;
        lossVal = exp(-r * curTime) * (1 - equivNList(curStep));
        totalLoss = totalLoss + lossVal * simulationDt;
        %accumLoss(curStep) = 100*totalLoss/curTime;
        accumLoss(curStep) = 100*totalLoss/finalTime;
    end

    ret = accumLoss;


end