function ret = costFunctionIntegral(x, beta, simulationDt)
    dataSize = size(x);
    finalStep = dataSize(2);
    ret = 0;
    for curStep = 1 : finalStep 
        curTime = (curStep - 1) * simulationDt;
        curCost = costFunction(x(:, curStep), beta(:, curStep), curTime);
        ret = ret + curCost * simulationDt; 
    end


end