addpath("..\Common\")
betaNList = [];
betaNList = [0.53];
loadFileName = ".\ContinuousSensitivityResultsWithEndo\beta_N_0.53_useArmijo.mat";
%loadFileName = "..\DiscreteTime\DiscreteSensitivityResultsWithEndo\beta_N_0.53_newT1.mat";
%loadFileName = ".\noEndoTry.mat";
%loadFileName = ".\ContinuousSensitivityResultsWithEndo\beta_N_0.53_new_with_descrete_guess.mat";
%loadFileName = "";
% for curBetaN = betaNList
%     paramMap = containers.Map(["betaN"], [curBetaN]);
%     saveFileName = ".\ContinuousSensitivityResultsWithEndo\beta_N_" + curBetaN + "_smoothed_discrete_guess.mat";
%     %saveFileName = "";
%     continuousOptimization(paramMap, saveFileName, loadFileName);
% end
% return;
loadFileName = "..\DiscreteTime\DiscreteSensitivityResultsWithEndo\beta_W_0.35.mat";
betaWList = [0.2,0.3,0.376, 0.4,0.5,0.6];
betaWList = [0.4];
for curBetaW = betaWList
    paramMap = containers.Map(["betaW"], [curBetaW]);
    %saveFileName = ".\ContinuousSensitivityResultsWithEndo\beta_W_" + curBetaW + "_smoothed_discrete_guess.mat";
    saveFileName = ".\ContinuousSensitivityResultsWithEndo\beta_W_" + curBetaW + "_NEW.mat";
    continuousOptimization(paramMap, saveFileName, "");
end
return

TvList = [500, 520, 540, 560, 580, 600];
TvList = [500, 600];
eulerConst = 0.577216;
sigmaTv = 44.74;
getMuTv = @(t)(t + eulerConst * sigmaTv);
for curTv = TvList 
    muTv = getMuTv(curTv);
    paramMap = containers.Map(["muTv"],[muTv]);
    saveFileName = "./ContinuousSensitivityResultsWithEndo/Tv_" + curTv + "new.mat";
    continuousOptimization(paramMap, saveFileName, "");
end



