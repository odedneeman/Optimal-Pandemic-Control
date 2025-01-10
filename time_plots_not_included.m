% Here we plot time charts, most not included in the paper.
% In each chart we have 5 plots: Employment, Effective R, Span of the disease,
% Infections per 100K (rate), and death rate per 100K. All as function of time.
%
% We have 4 kinds of charts- policy menu, stepwise, robustness and
% tracking.
% Each type (besides policy menu) can be made for 4 kinds of parameters-
% betaN, betaW, Tv, delta.
% 
% Policy menu- The disease has the benchmark parameters. We have 4 lines:
%   1. The stepwise optimal policy
%   2. The continuous optimal policy
%   3. No intervention
%   4. Full lockdown- that is employment at 75%
%
% Stepwise- the optimal stepwise policy, that's the sensitiviy analysis.
%
% Robustness- assuming the disease has the benchmark parameters and
% planning policy accordingly (i.e. stepwise optimal for benchmark), but
% one of the parameter has a different value.
%
% Tracking- effective R tracking.

clear;
close all;
addpath("Common\");

chart = 4; % 1 - policy menu. 2 - stepwise. 3 - robustness. 4- tracking.
valueType = 2; % 0- policy menu. 1- betaN. 2- betaW. 3- Tv. 4- delta
plotNameSuff = "2";
writeFile = 0;

% Define euler constant (will use it later):
eulerConst = 0.577216;
runNames = {};
betaNRunNames = {"DiscreteTime\DiscreteSensitivityResultsWithEndo\beta_N_0.53.mat"};
betaWRunNames = {"DiscreteTime\DiscreteSensitivityResultsWithEndo\beta_N_0.53.mat"};
tVRunNames = {"DiscreteTime\DiscreteSensitivityResultsWithEndo\beta_N_0.53.mat"};
deltaRunNames = {"DiscreteTime\DiscreteSensitivityResultsWithEndo\beta_N_0.53.mat"};
betaNLegend = {"Benchmark"};
betaWLegend = {"Benchmark"};
tVLegend = {"Benchmark"};
deltaLegend = {"Benchmark"};
%colors = {'k', 'g', 'r', 'b', 'b--', 'r--', 'g--', 'k--', 'g-.', 'r-.', 'b-.', 'k-.'};
betaNColors = {'k', 'r', [1, 0.64, 0], 'b'};
betaWColors = {'k', 'b', [0.93, 0.51, 0.93], [1, 0.64, 0], 'r'};
tVColors = {'k', 'b', 'r'};
deltaColors = {'k', 'b', [0.93, 0.51, 0.93], [1, 0.64, 0], 'r'};

betaNList = {0.4, 0.45, 0.6};
betaWList = {0.3, 0.35, 0.4, 0.45};
tvList = {500, 600};
deltaList = {0.004, 0.006, 0.01, 0.012};

if (chart == 1 && valueType ~= 0) || (chart ~= 1 && valueType == 0)
    disp("Error in valueType")
    return;
end

%POLICY MENU:
if chart == 1
    runNames{end+1} = "DiscreteTime\DiscreteSensitivityResultsWithEndo\beta_N_0.53.mat";
    runNames{end+1} = "ContinuousTime\ContinuousSensitivityResultsWithEndo\beta_N_0.53.mat";
    runNames{end+1} = 'ContinuousTime\ContinuousSensitivityResultsWithEndo\No_Intervention.mat';
    runNames{end+1} = 'ContinuousTime\ContinuousSensitivityResultsWithEndo\Full_Lockdown.mat';
    legendNames = {"Stepwise Optimal Policy", "Continuous Optimal Policy", "No Intervention", "Full Lockdown"};
    colorNames = {'k', 'g', 'r', 'b'};
    plotName = "policy_menu";
end


%STEPWISE:
if chart == 2
    plotName = "sensitivity";
    for i = 1:numel(betaNList)
        beta = betaNList{i};
        betaNRunNames{end+1} = "DiscreteTime\DiscreteSensitivityResultsWithEndo\beta_N_" + beta + ".mat";
        betaNLegend{end+1} = "betaN = " + beta;
    end
    for i = 1:numel(betaWList)
        beta = betaWList{i};
        if beta == 0.45
            %beta = "0.45_new4";
        end
        betaWRunNames{end+1} = "DiscreteTime\DiscreteSensitivityResultsWithEndo\beta_W_" + beta + ".mat";
        betaWLegend{end+1} = "betaW = " + beta;
    end
    for i = 1:numel(tvList)
        tvVal = tvList{i};
        tVRunNames{end+1} = "DiscreteTime\DiscreteSensitivityResultsWithEndo\Tv_" + tvVal + ".mat";
        tVLegend{end+1} = "Vaccine Arrival = " + tvVal;
    end
    for i = 1:numel(deltaList)
        deltaVal = deltaList{i};
        deltaRunNames{end+1} = "DiscreteTime\DiscreteSensitivityResultsWithEndo\delta_" + deltaVal + ".mat";
        deltaLegend{end+1} = "delta = " + deltaVal;
    end
end


%ROBUSTNESS
if chart == 3
    legendNames = {};
    plotName = "robustness";
    for i = 1:numel(betaNList)
        beta = betaNList{i};
        betaNRunNames{end+1} = "DiscreteTime\DiscreteRobustnessResults\beta_N_benchmark_0.53_reality_" + beta + ".mat";
        betaNLegend{end+1} = "betaN = " + beta;
    end
    for i = 1:numel(betaWList)
        beta = betaWList{i};
        betaWRunNames{end+1} = "DiscreteTime\DiscreteRobustnessResults\beta_W_benchmark_0.376_reality_" + beta + ".mat";
        betaWLegend{end+1} = "betaW = " + beta;
    end
    for i = 1:numel(tvList)
        tvVal = tvList{i};
        tVRunNames{end+1} = "DiscreteTime\DiscreteRobustnessResults\Tv_benchmark_540_reality_" + tvVal + ".mat";
        tVLegend{end+1} = "Vaccine Arrival = " + tvVal;
    end
end

%TRACKING:
if chart == 4
    legendNames = {};
    tvList = {};
    plotName = "feedback_control";
    for i = 1:numel(betaNList)
        beta = betaNList{i};
        betaNRunNames{end+1} = "FeedbackControl\matResults\betaNEffectiveR_" + num2str(beta) + "_New.mat";
        betaNLegend{end+1} = "betaN = " + beta;
    end
    for i = 1:numel(betaWList)
        beta = betaWList{i};
        betaWRunNames{end+1} = "FeedbackControl\matResults\betaWEffectiveR_" + num2str(beta) + "_New.mat";
        betaWLegend{end+1} = "betaW = " + beta;
    end
    
end



if valueType == 1
    plotName = plotName + "_betaN";
    runNames = betaNRunNames;
    legendNames = betaNLegend;
    colorNames = betaNColors;
end

if valueType == 2
    plotName = plotName + "_betaW";
    runNames = betaWRunNames;
    legendNames = betaWLegend;
    colorNames = betaWColors;
end

if valueType == 3
    plotName = plotName + "_Tv";
    runNames = tVRunNames;
    legendNames = tVLegend;
    colorNames = tVColors;
end

if valueType == 4
    plotName = plotName + "_delta";
    runNames = deltaRunNames;
    legendNames = deltaLegend;
    colorNames = deltaColors;
end

numRuns = numel(runNames);


figure(1);
hold on;
figure(2);
hold on;
figure(3);
hold on;
figure(4);
hold on;
figure(5);
hold on;
for run = 1:numRuns
    load(runNames{run});
    tV = round(muTv - eulerConst * sigmaTv);
    [betaListVac, xListVac, nListVac] = getVaccineState(xList, betaList, nList, tV, simulationDt, finalStep);
    figure(1);
    plot(timeList, nListVac, color=colorNames{run}); % Employment
    EffectiveR = (betaListVac / gamma) .* xListVac(1, :);
    EffectiveR(tV*100:end) = 0; % Forcing R-Effective to be 0 after vaccine arrival
    figure(2);
    plot(timeList(1:tV*100-1), EffectiveR(1:tV*100-1), color=colorNames{run}); % R-Effective until the vaccine
    figure(3);
    plot(timeList, 1 - xListVac(1, :), color=colorNames{run}); % 1-S
    figure(4);
    plot(timeList, xListVac(3, :) * 100000, color=colorNames{run}); % I per 100K
    dDot = delta * theta * xListVac(4, :);
    figure(5);
    plot(timeList, dDot * 100000, color=colorNames{run}); % dDot per 100k
end

figure(1)
title("Employment");
xlim([0,finalTime]);
ylim([0.68,1]);
legend(runNames);

figure(2)
title("Effective R");
xlim([0,finalTime]);
ylim([0.67, 1.2]);
plot(timeList, ones(size(timeList)), 'k:', 'LineWidth', 1.5); % Adding a line for 1
legend([runNames, "R-Effective=1"]);

figure(3)
title("Span of the Disease");
xlim([0,finalTime]);
ylim([0,0.36]);
legend(runNames);

figure(4)
title("Infections per 100K");
xlim([0,finalTime]);
ylim([1, 1000]);
%yticks(10.^(0:0.5:5));  
legend(runNames);
set(gca, 'YScale', 'log'); % Logarithmic scale
ytickformat('%,.0f');

figure(5)
title("Death Rate per 100K");
xlim([0,finalTime]);
ylim([10^-3, 1])
legend(runNames);
set(gca, 'YScale', 'log');
ytickformat('%,.0f');


% Create a new figure
combinedFig = figure('Position', [100, 100, 500, 125]);

% Define the layout of the combined plot
rows = 1;
cols = 5;

% Loop through each subplot
for i = 1:5
    % Extract the axes from the i-th figure
    ax = subplot(rows, cols, i);
    
    % Get the handle of the i-th figure
    figHandle = figure(i);
    
    % Get the axes from the i-th figure
    figAxes = get(figHandle, 'CurrentAxes');
    
    % Copy the children (plots, lines, etc.) from the i-th axes to the new axes
    copyobj(allchild(figAxes), ax);
    
    % Set the title of the subplot if needed
    titleString = get(get(figAxes, 'Title'), 'String');
    title(ax, titleString);
    
    % Copy the y-axis scale properties from the original axes to the new axes
    set(ax, 'YScale', get(figAxes, 'YScale'));
    
    % Set the x-axis limits to match the original plot
    xlim(ax, xlim(figAxes));
    ylim(ax, ylim(figAxes));
    set(ax, 'FontSize', 4);

    % Close the original figure to avoid cluttering the workspace
    close(figHandle);
end
legend(legendNames);
subplot(1, 5, 4);
yticks(10.^(0:3));

return
fig = gcf; % Get current figure handle
fig.Position(3:4) = [500 125]; % Set width and height in pixels
fileName = plotName + plotNameSuff + ".png";
%saveas(gca, fileName);

if writeFile == 0
    return
end

print(fileName, '-dpng', '-r600'); % 600 dpi resolution


% Read the PNG image
imageData = imread(fileName);

% Specify the region of interest (ROI) to crop
% Define the coordinates of the bounding box [xmin ymin width height]
roi = [300 0 2570 800];

% Crop the image using the ROI
croppedImage = imcrop(imageData, roi);

% Save the cropped image as a new PNG file
imwrite(croppedImage, fileName);
return

