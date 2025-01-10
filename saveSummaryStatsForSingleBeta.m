clear;
addpath("Common\");

paramType = "betaN";
runType = "stepwise";
% Define euler constant (will use it later):
eulerConst = 0.577216;

% The 5 summary statistics:
variableNames = {'Death Toll', 'GDP Loss', 'Span', 'Welfare loss', 'J Value'};
numVariables = numel(variableNames);
stats = zeros(1, numVariables);
betaWList = {0.3, 0.35, 0.376, 0.4, 0.45};
%betaWList = {0.374};
betaNList = {0.4, 0.45, 0.6};
betaNList = {0.41, 0.42, 0.42, 0.43, 0.44, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59};
betaWList = {0.32, 0.384};
betaNList = {0.55};
paramList = betaWList;
if paramType == "betaN"
    paramList = betaNList;
end


if(runType == "stepwise")
    for i = 1:numel(paramList)
        betaValue = paramList{i};
        if paramType == "betaW"
            load("DiscreteTime\DiscreteSensitivityResultsWithEndo\beta_W_" + betaValue);
        end
        if paramType == "betaN"
            load("DiscreteTime\DiscreteSensitivityResultsWithEndo\beta_N_" + betaValue);
        end
        % we need tV because we calculate the measures in time tV
        tV = round(muTv - eulerConst * sigmaTv);
        tvStep = round(tV / simulationDt);
        [betaListVac, xListVac, nListVac] = getVaccineState(xList, betaList, nList, tV, simulationDt, finalStep);
        stats(1) = xListVac(5, tvStep) * 100000; % Death toll
        stats(2) = 100 * getGDPLoss(xListVac, nListVac, simulationDt, finalTime); % GDP loss
        stats(3) = 100 * (1- xListVac(1, tvStep)); %  Span of the disease
        stats(4) = 100 * getEquivalentLoss(xListVac, betaListVac, simulationDt, finalStep) / finalTime; % Harm
        stats(5) = costFunctionIntegral(xList, betaList, simulationDt); % J
        % [equivLoss, integrand, equivNlist] = getEquivalentLossForTests(xListVac, betaListVac, simulationDt, finalStep);
        % stats(4) = (100 * equivLoss) / finalTime;
        % stats(6) = equivNlist;
        % stats(7) = integrand;
        disp("run " + num2str(i) + "");
        save("singleSummaryStats\" + paramType + "_" + betaValue + "_Stats.mat" ,"stats");
    end
end

if(runType == "feedback")
    for i = 1:numel(paramList)
        betaValue = paramList{i};
        if paramType == "betaW"
            load("FeedbackControl\matResults\betaWEffectiveR_" + betaValue + "_New");
        end
        if paramType == "betaN"
            load("FeedbackControl\matResults\betaNEffectiveR_" + betaValue + "_New");
        end
        % we need tV because we calculate the measures in time tV
        tV = round(muTv - eulerConst * sigmaTv);
        tvStep = round(tV / simulationDt);
        [betaListVac, xListVac, nListVac] = getVaccineState(xList, betaList, nList, tV, simulationDt, finalStep);
        stats(1) = xListVac(5, tvStep) * 100000; % Death toll
        stats(2) = 100 * getGDPLoss(xListVac, nListVac, simulationDt, finalTime); % GDP loss
        stats(3) = 100 * (1- xListVac(1, tvStep)); %  Span of the disease
        stats(4) = 100 * getEquivalentLoss(xListVac, betaListVac, simulationDt, finalStep) / finalTime; % Harm
        stats(5) = costFunctionIntegral(xList, betaList, simulationDt); % J
        % [equivLoss, integrand, equivNlist] = getEquivalentLossForTests(xListVac, betaListVac, simulationDt, finalStep);
        % stats(4) = (100 * equivLoss) / finalTime;
        % stats(6) = equivNlist;
        % stats(7) = integrand;
        disp("run " + num2str(i) + "");
        save("singleSummaryStats\feedbackControl\" + paramType + "_" + betaValue + "_Stats.mat", "stats");
    end
end

if(runType == "robustness")
    for i = 1:numel(paramList)
        betaValue = paramList{i};
        if paramType == "betaW"
            load("DiscreteTime\DiscreteRobustnessResults\beta_W_benchmark_0.376_reality_" + betaValue + ".mat");
        end
        if paramType == "betaN"
            load("DiscreteTime\DiscreteRobustnessResults\beta_N_benchmark_0.53_reality_" + betaValue + ".mat");
        end
        % we need tV because we calculate the measures in time tV
        tV = round(muTv - eulerConst * sigmaTv);
        tvStep = round(tV / simulationDt);
        [betaListVac, xListVac, nListVac] = getVaccineState(xList, betaList, nList, tV, simulationDt, finalStep);
        stats(1) = xListVac(5, tvStep) * 100000; % Death toll
        stats(2) = 100 * getGDPLoss(xListVac, nListVac, simulationDt, finalTime); % GDP loss
        stats(3) = 100 * (1- xListVac(1, tvStep)); %  Span of the disease
        stats(4) = 100 * getEquivalentLoss(xListVac, betaListVac, simulationDt, finalStep) / finalTime; % Harm
        stats(5) = costFunctionIntegral(xList, betaList, simulationDt); % J
        % [equivLoss, integrand, equivNlist] = getEquivalentLossForTests(xListVac, betaListVac, simulationDt, finalStep);
        % stats(4) = (100 * equivLoss) / finalTime;
        % stats(6) = equivNlist;
        % stats(7) = integrand;
        disp("run " + num2str(i) + "");
        save("singleSummaryStats\robustness\" + paramType + "_" + betaValue + "_Stats.mat", "stats");
    end
end

if(runType == "continuous")
    for i = 1:numel(paramList)
        betaValue = paramList{i};
        if paramType == "betaW"
            load("ContinuousTime\ContinuousSensitivityResultsWithEndo\beta_W_" + betaValue);
            if betaValue == 0.376
               load("ContinuousTime\ContinuousSensitivityResultsWithEndo\beta_N_0.53_smoothed_discrete_guess.mat");
            end
            if betaValue == 0.3
               load("ContinuousTime\ContinuousSensitivityResultsWithEndo\beta_W_0.3_smoothed_discrete_guess.mat");
            end
            if betaValue == 0.35
                load("ContinuousTime\ContinuousSensitivityResultsWithEndo\beta_W_0.35_NEW.mat");
            end
            if betaValue == 0.4
                load("ContinuousTime\ContinuousSensitivityResultsWithEndo\beta_W_0.4_NEW.mat");
            end
        end
        if paramType == "betaN"
            load("ContinuousTime\ContinuousSensitivityResultsWithEndo\beta_N_" + betaValue);
        end
        % we need tV because we calculate the measures in time tV
        tV = round(muTv - eulerConst * sigmaTv);
        tvStep = round(tV / simulationDt);
        [betaListVac, xListVac, nListVac] = getVaccineState(xList, betaList, nList, tV, simulationDt, finalStep);
        stats(1) = xListVac(5, tvStep) * 100000; % Death toll
        stats(2) = 100 * getGDPLoss(xListVac, nListVac, simulationDt, finalTime); % GDP loss
        stats(3) = 100 * (1- xListVac(1, tvStep)); %  Span of the disease
        stats(4) = 100 * getEquivalentLoss(xListVac, betaListVac, simulationDt, finalStep) / finalTime; % Harm
        stats(5) = costFunctionIntegral(xList, betaList, simulationDt); % J
        % [equivLoss, integrand, equivNlist] = getEquivalentLossForTests(xListVac, betaListVac, simulationDt, finalStep);
        % stats(4) = (100 * equivLoss) / finalTime;
        % stats(6) = equivNlist;
        % stats(7) = integrand;
        disp("run " + num2str(i) + "");
        save("singleSummaryStats\continuous\" + paramType + "_" + betaValue + "_StatsNew.mat", "stats");
    end
end
% paramList = {0.4, 0.45, 0.6};
% paramType = "betaN";
% for i = 1:numel(paramList)
%         betaValue = paramList{i};
%         load("singleSummaryStats\" + paramType + "_" + betaValue + "_Stats.mat");
%         save("singleSummaryStats\full_ones\" + paramType + "_" + betaValue + "_Stats.mat")
%         save("singleSummaryStats\" + paramType + "_" + betaValue + "_Stats.mat", "stats");
% end