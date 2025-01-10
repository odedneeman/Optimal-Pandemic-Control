% Figure 3.
% We plot here 3 plots. Each is a time series that corresponds to another
% betaW value. In each plot we present the stepwise, continuous and fixed
% policy that are optimal for this betaW value.


clear;
close all;
addpath("Common\");

fileName = "";
writeFile = 0;
lineWidth = 2;
getLegend = 1;
axisLabels = 1;

%contColor = [73, 159, 104]/255; % (Shamrock) Green
%fixedColor = [243, 201, 139]/255; % (Sunset) Orange
%stepwiseColor = [25, 123, 189]/255; % (Honolulu) Blue
betaW3Color = [47, 151, 193]/255;
betaW376Color = 'k';
betaW45Color  = 'r';

contStyle = '--';
stepwiseStyle = '-';
fixedStyle = '-.';
eulerConst = 0.577216;

paramValues = {0.3, 0.376, 0.45};

load("ContinuousTime\ContinuousSensitivityResultsWithEndo\beta_N_0.53.mat");

%tV = round(muTv - eulerConst * sigmaTv);
%[betaListVac, xListVac, nListVac] = getVaccineState(xList, betaList, nList, tV, simulationDt, finalStep);
nReduction = 100*(1 - nList);
figure(2);
plot(timeList, nReduction, 'Color', betaW376Color, 'LineWidth', lineWidth, 'LineStyle', contStyle);
hold on


load("DiscreteTime\DiscreteSensitivityResultsWithEndo\beta_N_0.53.mat");
nReduction = 100*(1 - nList);
figure(2);
plot(timeList, nReduction, 'Color', betaW376Color, 'LineWidth', lineWidth, 'LineStyle', stepwiseStyle);
xlim([0, 635]);
if axisLabels == 1
    xlabel('Time (days)');
    ylabel('Employment Reduction (%)');
end
ylim([0, 26]);
%title('betaW = 0.376');


load("ContinuousTime\statsPerN_Benchmark.mat");

nLevels = statsPerN(1, :);
harm = statsPerN(4, :);
[minHarm, minIndex] = min(harm);
fixedNHarm = minHarm;
bestFixedPol = 100*(1 - nLevels(minIndex));

figure(2)
plot(xlim, [bestFixedPol, bestFixedPol], 'Color', betaW376Color, 'LineWidth', lineWidth, 'LineStyle', fixedStyle);




% load("ContinuousTime\statsPerN_betaW_0.45.mat");
% nLevels = statsPerParam(1, :);
% harm = statsPerParam(4, :);
% [minHarm, minIndex] = min(harm);
% fixedNHarm = minHarm;
% bestFixedN = nLevels(minIndex);
%plot(xlim, [bestFixedN, bestFixedN]);






% load("FeedbackControl\matResults\betaWEffectiveR_0.45_New.mat");
% nReduction = 100*(1 - nList);
% figure(2)
% plot(timeList, nReduction, 'Color', 'g', 'LineWidth', lineWidth);

load("ContinuousTime\ContinuousSensitivityResultsWithEndo\beta_W_0.3.mat");

%tV = round(muTv - eulerConst * sigmaTv);
%[betaListVac, xListVac, nListVac] = getVaccineState(xList, betaList, nList, tV, simulationDt, finalStep);
nReduction = 100*(1 - nList);
figure(1);
plot(timeList, nReduction,  'Color', betaW3Color, 'LineWidth', lineWidth, 'LineStyle', contStyle);
hold on

load("DiscreteTime\DiscreteSensitivityResultsWithEndo\beta_W_0.3.mat");
nReduction = 100*(1 - nList);
figure(1);
plot(timeList, nReduction, 'Color', betaW3Color, 'LineWidth', lineWidth, 'LineStyle', stepwiseStyle);
xlim([0, 635]);
if axisLabels == 1
    xlabel('Time (days)');
    ylabel('Employment Reduction (%)');
end
ylim([0, 26]);
%title('betaW = 0.3');


load("ContinuousTime\statsPerN_betaW_0.3.mat");

nLevels = statsPerParam(1, :);
harm = statsPerParam(4, :);
[minHarm, minIndex] = min(harm);
fixedNHarm = minHarm;
bestFixedPol = 100*(1 - nLevels(minIndex));

figure(1)
plot(xlim, [bestFixedPol, bestFixedPol],  'Color', betaW3Color, 'LineWidth', lineWidth, 'LineStyle', fixedStyle);


load("ContinuousTime\ContinuousSensitivityResultsWithEndo\beta_W_0.45.mat");
nReduction = 100*(1 - nList);
figure(3)
plot(timeList, nReduction,  'Color', betaW45Color, 'LineWidth', lineWidth, 'LineStyle', contStyle);
hold on

load("DiscreteTime\DiscreteSensitivityResultsWithEndo\beta_W_0.45.mat");
nReduction = 100*(1 - nList);
figure(3);
plot(timeList, nReduction,  'Color', betaW45Color, 'LineWidth', lineWidth, 'LineStyle', stepwiseStyle);
xlim([0, 635]);
if axisLabels == 1
    xlabel('Time (days)');
    ylabel('Employment Reduction (%)');
end
ylim([0, 26]);
%title('betaW = 0.45');

load("ContinuousTime\statsPerN_betaW_0.45.mat");

nLevels = statsPerParam(1, :);
harm = statsPerParam(4, :);
[minHarm, minIndex] = min(harm);
fixedNHarm = minHarm;
bestFixedPol = 100*(1 - nLevels(minIndex));

figure(3)
plot(xlim, [bestFixedPol, bestFixedPol],  'Color', betaW45Color, 'LineWidth', lineWidth, 'LineStyle', fixedStyle);


combinedFig = figure('Position', [100, 100, 500, 125]);
rows = 1;
cols = 3;


% Loop through each subplot
for i = 1:3
    % Extract the axes from the i-th figure
    ax = subplot(rows, cols, 4-i);
    
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
    set(gca, 'FontSize', 12); % Changes font size for all text in the axes
end


fig = gcf; % Get current figure handle
fig.Position(3:4) = [1200 350]; % Set width and height in pixels
%fig.Position(3:4) = [300 1200];
if writeFile == 1
    print(fileName, '-dpng', '-r600'); % 600 dpi resolution
end



if getLegend == 1
    figure;

    % Create invisible plots with the desired styles
    plot(NaN, NaN,  'Color', 'k', 'LineWidth', lineWidth, 'LineStyle', fixedStyle); 
    hold on;
    plot(NaN, NaN,  'Color', 'k', 'LineWidth', lineWidth, 'LineStyle', contStyle); 
    plot(NaN, NaN, 'Color', 'k', 'LineWidth', lineWidth, 'LineStyle', stepwiseStyle); 
    
    legendNames = {'Optimal Fixed Policy', 'Optimal Continuous Policy', 'Optimal Stepwise Policy'};
    
    % Manually create the legend
    legend(legendNames);
    
    % Adjust the legend position or other properties if needed
    legend('Location', 'best');
end