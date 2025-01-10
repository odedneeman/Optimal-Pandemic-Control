load(".\ContinuousSensitivityResultsWithEndo\beta_N_0.53.mat")
wValueList = zeros(1, finalStep);
for curStep = 1 : finalStep
    curTime = (curStep - 1) * simulationDt;
    wValueList(curStep) = costFunction(xList(:,curStep), ...
                        betaList(curStep), curTime);
end
equivNList = zeros(1, finalStep);
options = optimset('Display','off');
for curStep = 1 : finalStep
    nEquation = @(n)(-log(n) - 1/5*(1 - n.^5) - wValueList(curStep));
    curNValue = fsolve(nEquation, 0.5, options);
    equivNList(curStep) = curNValue;
end

plot(0:0.01:635, equivNList)
hold on
plot(0:0.01:635,nList)
legend("equivalent employment rate")
legend(["equivalent employment rate", "true employment rate"])

measureTime = 540; %tV
measureFinalStep = measureTime / simulationDt;
totalLoss = 0;
for curStep = 1 : measureFinalStep + 1
    curTime = (curStep - 1)*simulationDt;
    lossVal = exp(-r * curTime) * (1 - equivNList(curStep));
    totalLoss = totalLoss + lossVal * simulationDt;
end
tvStep = measureFinalStep + 2; % The vaccine has arrived
disp(totalLoss);


