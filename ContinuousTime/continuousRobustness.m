benchmarkInputFile = ".\ContinuousSensitivityResultsWithEndo\beta_N_0.53.mat";
addpath("..\Common\")
load(benchmarkInputFile);
benchmarkCost = min(costList(costList>0));
benchmarkNList = nList;
benchmarkState = xList; 

fig1LegendList = ["Benchmark"];
figure(1); % Fig 1: cost list for different values 
hold on; 
stem(betaW, benchmarkCost);
fig2LegendList = ["Benchmark"];
figure(2); % Fig 2: death toll
hold on;
plot(timeList, xList(5,:) * 100000);
fig3LegendList = ["Benchmark"];
figure(3); % Fig 3: (real) employment rate
plot(timeList, nList);
hold on;

% Provide a list of parameters  
%betaNList = [0.4, 0.45, 0.5, 0.55, 0.6, 0.65];
betaWList = [0.2, 0.3, 0.4, 0.5, 0.6];
%tvList = 500:20:600;
betaListWithEndo = zeros(1, finalStep);
for newParameter = tvList
    betaW = newParameter;

    xList = zeros(6, finalStep);
    x0 = [0.9999; 0.00005; 0.00005; 0; 0; 0];
    xList(:,1) = x0;
    for curStep = 1 : finalStep
        curTime = (curStep - 1) * simulationDt;
        curX = xList(:, curStep);

        curN = benchmarkNList(curStep);
        curBeta = getBetaFromN(curN, curTime);
        betaRange = getBetaRange(curTime, curX);
        if curBeta < betaRange(1)
            curBeta = betaRange(1);
        end
        if curBeta > betaRange(2)
            curBeta = betaRange(2);
        end
        betaListWithEndo(curStep) = curBeta;
        dx = seirdDynamics(curX, curBeta);
        nextX = curX + dx * simulationDt;
        if curStep <= finalStep - 1
            xList(:, curStep + 1) = nextX;
        end
        nList(curStep) = getNFromBeta(curBeta, curTime);
    end
    newCost = costFunctionIntegral(xList, betaListWithEndo, simulationDt);
    newCost = real(newCost);
    
    fig1LegendList(end + 1) = "\beta_W = " + string(betaW);
    fig2LegendList(end + 1) = "\beta_W = " + string(betaW);
    fig3LegendList(end + 1) = "\beta_W = " + string(betaW);
    figure(1);
    stem(betaW, newCost);
    figure(2); 
    plot(timeList, xList(5,:) * 100000);
    figure(3);
    plot(timeList, nList);
end
figure(1)
xlabel("Values of \beta_W")
ylabel("Values of the objective function")
legend(fig1LegendList);
figure(2)
legend(fig2LegendList);
xlabel("Time (days)")
ylabel("Death toll (per 100,000)")
figure(3)
legend(fig3LegendList);
xlabel("Time (days)")
ylabel("Actual employment rate")