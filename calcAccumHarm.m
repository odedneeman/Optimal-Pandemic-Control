% In the summary stats we calculate the accumulated harm only at the end.
% Here we do it for every step.

betaWList1 = {0.3, 0.4, 0.45};
betaNList1 = {0.4, 0.45, 0.6};
eulerConst = 0.577216;


for i=1:numel(betaWList1)
    paramValue = betaWList1{i};
    %load("DiscreteTime\DiscreteSensitivityResultsWithEndo\beta_W_" + paramValue + ".mat");
    %load("DiscreteTime\DiscreteRobustnessResults\beta_W_benchmark_0.376_reality_" + paramValue + ".mat");
    load("FeedbackControl\matResults\betaWEffectiveR_" + paramValue + "_New.mat");
    tV = round(muTv - eulerConst * sigmaTv);
    tvStep = round(tV / simulationDt);
    [betaListVac, xListVac, nListVac] = getVaccineState(xList, betaList, nList, tV, simulationDt, finalStep);
    accumLoss = getEquivalentLossAllSteps(xListVac, betaListVac, simulationDt, finalStep);
    save("accumHarm\tracking\betaW" + "_" + paramValue + ".mat", "accumLoss")
    disp("betaW=" + paramValue + " done");
end

for i=1:numel(betaNList1)
    paramValue = betaNList1{i};
    %load("DiscreteTime\DiscreteSensitivityResultsWithEndo\beta_N_" + paramValue + ".mat");
    %load("DiscreteTime\DiscreteRobustnessResults\beta_N_benchmark_0.53_reality_" + paramValue + ".mat");
    load("FeedbackControl\matResults\betaNEffectiveR_" + paramValue + "_New.mat");
    tV = round(muTv - eulerConst * sigmaTv);
    tvStep = round(tV / simulationDt);
    [betaListVac, xListVac, nListVac] = getVaccineState(xList, betaList, nList, tV, simulationDt, finalStep);
    accumLoss = getEquivalentLossAllSteps(xListVac, betaListVac, simulationDt, finalStep);
    save("accumHarm\tracking\betaN" + "_" + paramValue + ".mat", "accumLoss")
    disp("betaN=" + paramValue + " done");
end
