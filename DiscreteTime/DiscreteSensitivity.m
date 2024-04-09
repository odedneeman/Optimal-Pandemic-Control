addpath("..\Common\")
betaNList = [0.4, 0.5, 0.53, 0.55, 0.6, 0.65];
betaNList = [];
for curBetaN = betaNList
    paramMap = containers.Map(["betaN"], [curBetaN]);
    saveFileName = "./DiscreteSensitivityResultsWithEndo/beta_N_" + curBetaN + ".mat";
    discreteOptimization(paramMap, saveFileName);
end

betaWList = [0.2,0.3,0.376, 0.4,0.5,0.6];
betaWList = [0.376, 0.4,0.5,0.6];
for curBetaW = betaWList
    paramMap = containers.Map(["betaW"], [curBetaW]);
    saveFileName = "./DiscreteSensitivityResultsWithEndo/beta_W_" + curBetaW + ".mat";
    discreteOptimization(paramMap, saveFileName);
end
return;
TvList = [500, 520, 540, 560, 580, 600];
TvList = [];
eulerConst = 0.577216;
sigmaTv = 44.74;
getMuTv = @(t)(t + eulerConst * sigmaTv);
for curTv = TvList 
    muTv = getMuTv(curTv);
    paramMap = containers.Map(["muTv"],[muTv]);
    saveFileName = "./DiscreteSensitivityResultsWithEndo/Tv_" + curTv + ".mat";
    discreteOptimization(paramMap, saveFileName);
end