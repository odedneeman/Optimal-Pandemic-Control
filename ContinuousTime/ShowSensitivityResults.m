paramList = [0.4, 0.45, 0.5, 0.53, 0.55, 0.6, 0.65];
targetParam = "betaN";
fileNameAppendix = ".\ContinuousSensitivityResultsWithEndo\beta_N_";

figure(1);
hold on;
figure(2);
hold on;
figure(3);
hold on;

for curBetaNIndex = 1 : length(paramList)
    curParamValue = paramList(curBetaNIndex);
    loadFileName = fileNameAppendix + curParamValue + ".mat";
    load(loadFileName)
    costList = real(costList);
    finalCost = min(costList(costList > 0));
    figure(1) 
    stem(curParamValue, finalCost);
    figure(2)
    tList = 0:0.01:finalTime;
    infected = xList(3,:);
    plot(tList, infected * 100000);
    
    figure(3)
    death = xList(5,:);
    plot(tList, death * 100000)
    
end
legendList = zeros(1, length(paramList));
legendList = string(legendList);
for curIndex = 1 : length(paramList)
    legendList(curIndex) = targetParam + "=" + string(paramList(curIndex));
end

for figureIndex = 1:3
    figure(figureIndex)
    legend(legendList)
end
figure(1);
xlabel("value of " + targetParam);
ylabel("(locally) optimal cost function value")
figure(2);
xlabel("time (days)")
ylabel("amount of infected people per 100,000")
figure(3);
xlabel("time (days)")
ylabel("cumulative death per 100,000")