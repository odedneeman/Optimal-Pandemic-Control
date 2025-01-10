addpath("..\Common\")
% betaNList = [0.4, 0.45, 0.5, 0.53, 0.55, 0.6];
betaNList = [0.55];

for curBetaN = betaNList
    paramMap = containers.Map(["betaN"], [curBetaN]);
    saveFileName = "./DiscreteSensitivityResultsWithEndo/beta_N_" + curBetaN + "_NEW2.mat";
    discreteOptimization(paramMap, saveFileName);
    disp("betaN=" + curBetaN + " done")
end
% betaWList = [0.2, 0.3, 0.376, 0.4, 0.5, 0.6];
% betaWList = [0.376, 0.4, 0.5 ,0.6];
% betaWList = [0.384]; %new ones to add
% for curBetaW = betaWList
%     paramMap = containers.Map(["betaW"], [curBetaW]);
%     saveFileName = "./DiscreteSensitivityResultsWithEndo/beta_W_" + curBetaW + "_new3.mat";
%     discreteOptimization(paramMap, saveFileName);
%     disp("betaW=" + curBetaW + " done")
% end
% TvList = [500, 520, 540, 560, 600];
% TvList = [];
% eulerConst = 0.5772   16;
% sigmaTv = 44.74;
% getMuTv = @(t)(t + eulerConst * sigmaTv);
% for curTv = TvList 
%     muTv = getMuTv(curTv);
%     paramMap = containers.Map(["muTv"],[muTv]);
%     saveFileName = "./DiscreteSensitivityResultsWithEndo/Tv_" + curTv + "new.mat";
%     discreteOptimization(paramMap, saveFileName);
%     disp("tV=" + curTv + " done")
% end
% deltaList = [0.012];
% for curDelta = deltaList
%     paramMap = containers.Map(["delta"],[curDelta]);
%     saveFileName = "./DiscreteSensitivityResultsWithEndo/delta_" + curDelta + "_3.mat";
%     discreteOptimization(paramMap, saveFileName);
%     disp("delta=" + curDelta + " done")
% end

