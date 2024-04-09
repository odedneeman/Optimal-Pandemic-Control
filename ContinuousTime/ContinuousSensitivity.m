addpath("..\Common\")
betaNList = [];
%betaNList = [0.53];
loadFileName = ".\ContinuousSensitivityResultsWithEndo\beta_N_0.45_best.mat";
%loadFileName = ".\ContinuousSensitivityResultsWithEndo\beta_N_0.53.mat";
%loadFileName = ".\ContinuousSensitivityResultsWithEndo\beta_N_0.45.mat";
%loadFileName = "";
for curBetaN = betaNList
    paramMap = containers.Map(["betaN"], [curBetaN]);
    saveFileName = "./ContinuousSensitivityResultsWithEndo/beta_N_" + curBetaN + "_from_0.45.mat";
    %saveFileName = "";
    continuousOptimization(paramMap, saveFileName, loadFileName);
end

betaWList = [0.2,0.3,0.376, 0.4,0.5,0.6];
%betaWList = [];
for curBetaW = betaWList
    paramMap = containers.Map(["betaW"], [curBetaW]);
    saveFileName = "./ContinuousSensitivityResultsWithEndo/beta_W_" + curBetaW + ".mat";
    continuousOptimization(paramMap, saveFileName, "");
end

TvList = [500, 520, 540, 560, 580, 600];
%TvList = [];
eulerConst = 0.577216;
sigmaTv = 44.74;
getMuTv = @(t)(t + eulerConst * sigmaTv);
for curTv = TvList 
    muTv = getMuTv(curTv);
    paramMap = containers.Map(["muTv"],[muTv]);
    saveFileName = "./ContinuousSensitivityResultsWithEndo/Tv_" + curTv + ".mat";
    continuousOptimization(paramMap, saveFileName, "");
end

