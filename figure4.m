% This is figure 4 - which demonstrates the robustness of the feedback control
% This script also creates Figure S4, by switching paramType to "betaN"
close all;
robColor = [108, 0, 122] / 255; % purple
feedbackColor = [99, 193, 50] / 255; % green
markerColor = [99/255, 193/255, 50/255, 0.3];
getLegend = 1;
notCombined = 0;
extraFigure = 1;
fileName = "difference_from_opt_plot_betaW_31.7.png";
writeFile = 0;
paramType = "betaN";


benchmarkParam = 0.376;
if paramType == "betaN"
    benchmarkParam = 0.53;
end
markerSize = 30;
lineWidth = 1.5;

xlabelStr = "Actual baseline transmission rate, \beta_W";
if paramType == "betaN"
    xlabelStr = "Actual employment impact level, \beta_N";
end
betaWList = [0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.366, 0.37, 0.374, 0.376, 0.38, 0.384, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45];
betaNList = 0.4:0.01:0.6;
paramList = betaWList;
if paramType == "betaN"
    paramList = betaNList;
end

numOfRes = numel(paramList);
deathTollList = zeros(1, numOfRes);
GDPLossList = zeros(1, numOfRes);
harmList = zeros(1, numOfRes);

for i = 1:numOfRes
    param = paramList(i);
    load("singleSummaryStats\" + paramType + "_" + param + "_Stats.mat", "stats");
    deathToll = stats(1);
    GDPLoss = stats(2);
    span = stats(3);
    harm = stats(4);
    deathTollList(i) = deathToll;
    GDPLossList(i) = GDPLoss;
    harmList(i) = harm;
end

% Robustness:
    
load("ContinuousTime\statsPer" + paramType + "_Benchmark.mat");

betas = statsPerParam(1, :);
deathToll = statsPerParam(2, :);
GDPLoss = statsPerParam(3, :);
harm = statsPerParam(4, :);
span = statsPerParam(5, :);

deathTollDiffRobustness = zeros(1, numOfRes);
GDPLossDiffRobustness = zeros(1, numOfRes);
harmDiffRobustness = zeros(1, numOfRes);

for i = 1:numOfRes
    param = paramList(i);
    tolerance = 1e-6;  % Small tolerance value
    j = find(abs(betas - param) < tolerance);
    deathTollDiffRobustness(i) = deathToll(j) - deathTollList(i);
    GDPLossDiffRobustness(i) = GDPLoss(j) - GDPLossList(i);
    harmDiffRobustness(i) = harm(j) - harmList(i);
end

figure(1)
%plot(param_interp, deathTollDiffRobustness, 'LineWidth', 1.5, 'Color', robColor);
%scatter(betaWList, deathTollDiffRobustness, 'Color', robColor);
scatter(paramList, deathTollDiffRobustness, markerSize, 'filled', 'MarkerFaceColor', robColor, 'MarkerFaceAlpha', 0.7);
hold on

figure(2)
%plot(param_interp, GDPLossDiffRobustness, 'LineWidth', 1.5, 'Color', robColor);
scatter(paramList, GDPLossDiffRobustness, markerSize, 'filled', 'MarkerFaceColor', robColor, 'MarkerFaceAlpha', 0.7);
hold on

load("FeedbackControl\matResults\statsPer" + paramType + "EffectiveRControl.mat");

betas = statsPerParam(1, :);
deathToll = statsPerParam(2, :);
GDPLoss = statsPerParam(3, :);
harmFeedback = statsPerParam(4, :);
span = statsPerParam(5, :);


deathTollDiffTrack = zeros(1, numOfRes);
GDPLossDiffTrack = zeros(1, numOfRes);
harmDiffTrack = zeros(1, numOfRes);

for i = 1:numOfRes
    param = paramList(i);
    tolerance = 1e-6;  % Small tolerance value
    j = find(abs(betas - param) < tolerance);
    deathTollDiffTrack(i) = deathToll(j) - deathTollList(i);
    GDPLossDiffTrack(i) = GDPLoss(j) - GDPLossList(i);
    harmDiffTrack(i) = harmFeedback(j) - harmList(i);
end

figure(1)
%plot(param_interp, deathTollDiffTrack, 'LineWidth', 1.5, 'Color', feedbackColor);
scatter(paramList, deathTollDiffTrack, markerSize, 'filled', 'MarkerFaceColor', feedbackColor, 'MarkerFaceAlpha', 0.7);

figure(2)
%plot(param_interp, GDPLossDiffTrack, 'LineWidth', 1.5, 'Color', feedbackColor);
scatter(paramList, GDPLossDiffTrack, markerSize, 'filled', 'MarkerFaceColor', feedbackColor, 'MarkerFaceAlpha', 0.7);

figure(1)
%plot([0, 0], ylim, 'k--'); % Vertical line at x = 0
plot(xlim, [0, 0], 'k--'); % Horizontal line at y = 1
%ylabel("Death Toll - Difference from Optimal Policy (Death Toll in reality - Optimal)")
ylabel({"Difference in death toll vs", "optimum (in deaths per 100K)"});
xlabel(xlabelStr)
xlim([paramList(1),paramList(end)]);
ylim([-100, 200]);
if paramType == "betaN"
    ylim([-20, 150]);
end
plot([benchmarkParam, benchmarkParam] , ylim, 'b--'); % Vertical line at x = BenchmarkParam

figure(2)
%plot([0, 0], ylim, 'k--'); % Vertical line at x = 0
plot(xlim, [0, 0], 'k--'); % Horizontal line at y = 0
%ylabel("GDP Loss Difference from Optimal Policy, in % of annual GDP (GDP Loss in reality - Optimal)")
ylabel({"Difference in GDP loss vs optimum", "(% of annual GDP)"});
xlabel(xlabelStr)
xlim([paramList(1),paramList(end)]);
ylim([-9, 16]);
if paramType == "betaN"
    ylim([-8, 5]);
end
if notCombined == 0
    plot([benchmarkParam, benchmarkParam] , ylim, 'b--'); % Vertical line at x = BenchmarkParam
end

figure(3)
%plot(param_interp, harmDiffRobustness, 'LineWidth', 1.5, 'Color', robColor);
scatter(paramList, harmDiffRobustness, markerSize, 'filled', 'MarkerFaceColor', robColor, 'MarkerFaceAlpha', 0.7);
hold on
%plot(param_interp, harmDiffTrack, 'LineWidth', 1.5, 'Color', feedbackColor);
scatter(paramList, harmDiffTrack, markerSize, 'filled', 'MarkerFaceColor', feedbackColor, 'MarkerFaceAlpha', 0.7);
plot(xlim, [0, 0], 'k--'); % Horizontal line at y = 0
%legend("Optimal Benchmark Policy", "Feedback Control")
xlabel(xlabelStr)
%ylabel("Harm Difference from Optimal Policy, in % of equivalent consumption loss (Harm in reality - Optimal)")
ylabel({"Difference in welfare loss vs optimum", "(% of equivalent consumption loss)"});
xlim([paramList(1),paramList(end)]);
if paramType == "betaN"
    ylim([-0.5, 3.2]);
end
plot([benchmarkParam, benchmarkParam] , ylim, 'b--'); % Vertical line at x = BenchmarkParam

for i = 1:3

    figure(i)
    % Get the current x-axis ticks and labels
    xticks = get(gca, 'XTick');
    xticklabels = get(gca, 'XTickLabel');
    
    % Add a new tick
    new_tick = benchmarkParam;
    new_xticks = [xticks, new_tick];
    new_xticks = sort(new_xticks);

    % Add a new label for the new tick
    %new_label = '';
    %new_xticklabels = [xticklabels; {new_label}];
    
    % Set the new x-axis ticks and labels
    set(gca, 'XTick', new_xticks);
    %set(gca, 'XTickLabel', new_xticklabels);
end

if notCombined == 1
    return
end

combinedFig = figure('Position', [100, 100, 500, 125]);
% Define the layout of the combined plot
rows = 1;
cols = 3;


% Loop through each subplot
for i = 1:3
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
    %set(ax, 'FontSize', 4);
    
    % Copy the x-axis label from the original axes to the new axes
    xlabelString = get(get(figAxes, 'XLabel'), 'String');
    xlabel(ax, xlabelString);
    
    % Copy the y-axis label from the original axes to the new axes
    ylabelString = get(get(figAxes, 'YLabel'), 'String');
    ylabel(ax, ylabelString);
    set(ax, 'XTick', get(figAxes, 'XTick'));
    set(ax, 'XTickLabel', get(figAxes, 'XTickLabel'));
    
    % Copy the x-axis ticks and labels from the original axes to the new axes
    set(ax, 'XTick', get(figAxes, 'XTick'));
    set(ax, 'XTickLabel', get(figAxes, 'XTickLabel'));

    % Copy the y-axis ticks and labels from the original axes to the new axes
    set(ax, 'YTick', get(figAxes, 'YTick'));
    set(ax, 'YTickLabel', get(figAxes, 'YTickLabel'));

    % Copy the legend if it exists in the original figure
    leg = get(figAxes, 'Legend');
    if ~isempty(leg)
        legend(ax, get(leg, 'String'), 'Location', 'best', 'FontSize', 7);
    end
    set(ax, 'FontSize', 8);
    

    % Close the original figure to avoid cluttering the workspace
    close(figHandle);
    set(gca, 'FontSize', 10); % Changes font size for all text in the axes
end


%legend("Robustness (Benchmark Stepwise Policy)", "Effective R Tracking");
fig = gcf; % Get current figure handle
fig.Position(3:4) = [1000 300]; % Set width and height in pixels
if writeFile == 1
    print(fileName, '-dpng', '-r600'); % 600 dpi resolution
end


if extraFigure == 1
    figure;
    betaW3Color = [47, 151, 193]/255;
    betaW376Color = 'k';
    betaW45Color  = 'r';
    load("ContinuousTime\ContinuousSensitivityResultsWithEndo\beta_W_0.3.mat");
    nReduction = 100*(1 - nList);
    plot(timeList, nReduction,  'Color', betaW3Color, 'LineWidth', 1.5, 'LineStyle', '--');
    hold on
    load("ContinuousTime\ContinuousSensitivityResultsWithEndo\beta_N_0.53.mat");
    nReduction = 100*(1 - nList);
    plot(timeList, nReduction, 'Color', betaW376Color, 'LineWidth', 1.5, 'LineStyle', '--');
    load("ContinuousTime\ContinuousSensitivityResultsWithEndo\beta_W_0.45.mat");
    nReduction = 100*(1 - nList);
    plot(timeList, nReduction,  'Color', betaW45Color, 'LineWidth', 1.5, 'LineStyle', '--');
    load("FeedbackControl\matResults\betaWEffectiveR_0.3_New.mat")
    nReduction = 100*(1 - nList);
    %plot(timeList, nReduction, 'Color', [0, 1, 0, 0.4], 'LineWidth', 4);
    %plot(timeList, nReduction, 'Color', markerColor, 'LineWidth', 6);
    plot(timeList, nReduction, 'Color', betaW3Color, 'LineWidth', lineWidth);
    load("FeedbackControl\matResults\betaWEffectiveR_0.376_New.mat")
    nReduction = 100*(1 - nList);
    %plot(timeList, nReduction, 'Color', markerColor, 'LineWidth', 6);
    plot(timeList, nReduction, 'Color', betaW376Color, 'LineWidth', lineWidth);
    load("FeedbackControl\matResults\betaWEffectiveR_0.45_New.mat")
    nReduction = 100*(1 - nList);
    %plot(timeList, nReduction, 'Color', markerColor, 'LineWidth', 6);
    plot(timeList, nReduction, 'Color', betaW45Color, 'LineWidth', lineWidth);
    xlim([0,635]);
    xlabel("Time (days)");
    ylabel("Employment Reduction (%)");
    set(gca, 'FontSize', 10);
    box off;
    fig = gcf;
    fig.Position(3:4) = [400 250];

    if getLegend == 1
        figure;

        % Create invisible plots with the desired styles
        plot(NaN, NaN,  'Color', betaW3Color, 'LineWidth', 2); 
        hold on;
        plot(NaN, NaN,  'Color', betaW376Color, 'LineWidth', 2); 
        plot(NaN, NaN, 'Color', betaW45Color, 'LineWidth', 2); 
        plot(nan, nan, 'w');
        %hold on;
        plot(nan, nan, 'w');
        %plot(nan, nan, 'w');

        legendNames = {'Low  Basic  Disease  Transimission', 'Medium  Basic  Disease  Transimission', 'High  Basic  Disease  Transmission', 'Dashed  Line  -  Optimal  Continuous  Policy  (with true knowledge about the disease)', 'Continuous  Line  -  Policy Using  Feedback  Control  Algorithm  (without prior knowledge)'};        
        % Manually create the legend

        legend(legendNames);
        fig = gcf;
        fig.Position(3:4) = [700 250];
        figure;
        scatter(NaN, NaN, markerSize, 'filled', 'MarkerFaceColor', robColor, 'MarkerFaceAlpha', 0.7);
        hold on;
        scatter(NaN, NaN, markerSize, 'filled', 'MarkerFaceColor', feedbackColor, 'MarkerFaceAlpha', 0.7);
        legendNames1 = {"Optimal policy for baseline transmission rate (\beta_W) of 0.376", "Policy updated using feedback control algorithm"};
        if paramType == "betaN"
            legendNames1 = {"Optimal policy for employment impact level (\beta_N) of 0.53", "Policy updated using feedback control algorithm"};
        end
        legend(legendNames1);

    end
end