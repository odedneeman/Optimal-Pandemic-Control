clear;
benchmarkInputFile = ".\DiscreteSensitivityResultsWithEndo\beta_N_0.53.mat";
addpath("..\Common\")
load(benchmarkInputFile);
benchmarkTv = 540;
benchmarkCost = min(costList(costList>0));
benchmarkNList = getContinuousEmployment(t1, n1, finalTime, simulationDt);
benchmarkState = xList; 

%no plotting for now

% fig1LegendList = ["Benchmark"];
% figure(1); % Fig 1: cost list for different values 
% hold on; 
% stem(betaN, benchmarkCost);
% %stem(betaW, benchmarkCost);
% fig2LegendList = ["Benchmark"];
% figure(2); % Fig 2: death toll
% hold on;
% plot(timeList, xList(5,:) * 100000);
% fig3LegendList = ["Benchmark"];
% figure(3); % Fig 3: (real) employment rate
% plot(timeList, nList);
% hold on;

% Provide a list of parameters 
tvList = [500, 600];

for tV = tvList
    [betaListVac, xListVac, nListVac] = getVaccineState(xList, betaList, nList, tV, simulationDt, finalStep);
    % To save the workspace we want the results to be in xList, betaList,
    % nList, so we will backup the original ones and then write back:
    xListBackUp = xList;
    nListBackUp = nList;
    betaListBackUp = betaList;
    xList = xListVac;
    betaList = betaListVac;
    nList = nListVac;
    % We'll need euler constant and sigmaTv to calculate muTv with tV to save it
    % to the workspace as well:
    eulerConst = 0.577216;
    sigmaTv = 44.74;
    muTv = tV + eulerConst * sigmaTv;
    fileName = "./DiscreteRobustnessResults/Tv_benchmark_" + benchmarkTv + "_reality_" + tV + ".mat";
    save(fileName);
    xList = xListBackUp;
    nList = nListBackUp;
    betaList = betaListBackUp;
end
