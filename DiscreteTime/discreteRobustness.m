clear;
benchmarkInputFile = ".\DiscreteSensitivityResultsWithEndo\beta_N_0.53.mat";
%benchmarkInputFile = ".\DiscreteSensitivityResultsWithEndo\beta_W_0.376.mat";
addpath("..\Common\")
load(benchmarkInputFile);
%benchmarkBeta = 0.53;
%benchmarkBeta = 0.376;
benchmarkCost = min(costList(costList>0));
benchmarkNList = getContinuousEmployment(t1, n1, finalTime, simulationDt);
benchmarkState = xList; 

fig1LegendList = ["Benchmark"];
figure(1); % Fig 1: cost list for different values 
hold on; 
%stem(betaN, benchmarkCost);
%stem(betaW, benchmarkCost);
fig2LegendList = ["Benchmark"];
figure(2); % Fig 2: death toll
hold on;
plot(timeList, xList(5,:) * 100000);
fig3LegendList = ["Benchmark"];
figure(3); % Fig 3: (real) employment rate
plot(timeList, nList);
hold on;

% Provide a list of parameters  
betaNList = [0.4, 0.45, 0.5, 0.6];
%betaWList = [0.3, 0.35, 0.4, 0.45];
deltaList = [0.004, 0.006, 0.008, 0.01, 0.012];
parameterList = deltaList;
betaListWithEndo = zeros(1, finalStep);
for newParameter = parameterList
    %betaN = newParameter;
    %betaW = newParameter;
    delta = newParameter;
    xList = zeros(6, finalStep + 1);
    betaListWithEndo = zeros(1, finalStep + 1);
    x0 = [0.9999; 0.00005; 0.00005; 0; 0; 0];
    xList(:,1) = x0;
    for curStep = 1 : finalStep + 1
        curTime = (curStep - 1) * simulationDt;
        curX = xList(:, curStep);

        curN = benchmarkNList(curStep);
        curBeta = getBetaFromN(curN, curTime); % Here "hides" the use of the new betaN/W
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
        if curStep <= finalStep
            xList(:, curStep + 1) = nextX;
        end
        nList(curStep) = getNFromBeta(curBeta, curTime);
    end
    newCost = costFunctionIntegral(xList, betaListWithEndo, simulationDt);
    newCost = real(newCost);
    betaList = betaListWithEndo;
    %fileName = "./DiscreteRobustnessResults/beta_N_benchmark_" + benchmarkBeta + "_reality_" + betaN + ".mat";
    %fileName = "./DiscreteRobustnessResults/beta_W_benchmark_" + benchmarkBeta + "_reality_" + betaW + ".mat";
    % fileName = "./DiscreteRobustnessResults/delta_benchmark_" + 0.008 + "_reality_" + delta + ".mat";
    % save(fileName);
    % fig1LegendList(end + 1) = "\beta_N = " + string(betaN);
    % fig2LegendList(end + 1) = "\beta_N = " + string(betaN);
    % fig3LegendList(end + 1) = "\beta_N = " + string(betaN);
    % % fig1LegendList(end + 1) = "\beta_W = " + string(betaW);
    % % fig2LegendList(end + 1) = "\beta_W = " + string(betaW);
    % % fig3LegendList(end + 1) = "\beta_W = " + string(betaW);
    % figure(1);
    % stem(betaN, newCost);
    % %stem(betaW, newCost);
    % figure(2); 
    % plot(timeList, xList(5,:) * 100000);
    % figure(3);
    % plot(timeList, nList);

end
% figure(1)
% xlabel("Values of \beta_N")
% %xlabel("Values of \beta_W")
% ylabel("Values of the objective function")
% legend(fig1LegendList);
% figure(2)
% legend(fig2LegendList);
% xlabel("Time (days)")
% ylabel("Death toll (per 100,000)")
% figure(3)
% legend(fig3LegendList);
% xlabel("Time (days)")
% ylabel("Actual employment rate")