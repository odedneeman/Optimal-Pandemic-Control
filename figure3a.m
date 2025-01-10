

% This is the top panel of Figure 3 
% The fillings were added manually.

close all;
load("ContinuousTime\statsPerbetaN_Benchmark.mat")

betaW3Color = [47, 151, 193]/255;
betaW376Color = 'k';
betaW45Color  = 'r';

getLegend = 1;
lineWidth = 2;

runTypes = {'No Policy', 'Optimal Fixed Policy', 'Optimal Stepwise Policy', 'Optimal Continuous Policy', 'Maximum Employment Reduction'};
fullLockdownIndex = 15;
noInterventionIndex = 65;

writeXTicks = 0;


harm_376_full_lockdown = real(statsPerN(4,fullLockdownIndex));
harm_376_no_intervention = real(statsPerN(4,noInterventionIndex));

load("ContinuousTime\statsPerN_betaW_0.45.mat");
harm_45_full_lockdown = statsPerParam(4, fullLockdownIndex);
harm_45_no_intervention = statsPerParam(4, noInterventionIndex);

load("ContinuousTime\statsPerN_betaW_0.3.mat");
harm_3_full_lockdown = statsPerParam(4, fullLockdownIndex);
harm_3_no_intervention = statsPerParam(4, noInterventionIndex);

load("ContinuousTime\statsPerN_Benchmark.mat");

nLevels = statsPerN(1, :);
harm = statsPerN(4, :);
[minHarm, minIndex] = min(harm);
harm_376_fixed = minHarm;
%bestFixedPol = 100*(1 - nLevels(minIndex));


load("ContinuousTime\statsPerN_betaW_0.3.mat");

nLevels = statsPerParam(1, :);
harm = statsPerParam(4, :);
[minHarm, minIndex] = min(harm);
harm_3_fixed = minHarm;
%bestFixedPol = 100*(1 - nLevels(minIndex));

load("ContinuousTime\statsPerN_betaW_0.45.mat");

nLevels = statsPerParam(1, :);
harm = statsPerParam(4, :);
[minHarm, minIndex] = min(harm);
harm_45_fixed = minHarm;
%bestFixedPol = 100*(1 - nLevels(minIndex));

load("singleSummaryStats\betaW_0.376_Stats.mat");
harm_376_stepwise = stats(4);

load("singleSummaryStats\betaW_0.3_Stats.mat");
harm_3_stepwise = stats(4);

load("singleSummaryStats\betaW_0.45_Stats.mat");
harm_45_stepwise = stats(4);

load("singleSummaryStats\continuous\betaW_0.376_Stats.mat");
harm_376_cont = stats(4);

load("singleSummaryStats\continuous\betaW_0.3_Stats.mat");
harm_3_cont = stats(4);

load("singleSummaryStats\continuous\betaW_0.45_Stats.mat");
harm_45_cont = stats(4);

betaW45Harms = [harm_45_no_intervention, harm_45_fixed, harm_45_stepwise, harm_45_cont, harm_45_full_lockdown];
betaW3Harms = [harm_3_no_intervention, harm_3_fixed, harm_3_stepwise, harm_3_cont, harm_3_full_lockdown];
betaW376Harms = [harm_376_no_intervention, harm_376_fixed, harm_376_stepwise, harm_376_cont, harm_376_full_lockdown];


% figure;
% 
% b45 = bar(betaW45Harms, LineWidth=lineWidth);
% b45.EdgeColor = betaW45Color;
% b45.FaceColor = 'none';
% hold on;
% bar(5, betaW45Harms(5), 'FaceColor', betaW45Color, 'EdgeColor', betaW45Color)

figure;

% Create the bar plot
b45 = bar(betaW45Harms, 'LineWidth', lineWidth);
b45.FaceColor = 'none';  % Make all bars empty
b45.EdgeColor = 'none';  % Remove the default edges
hold on;

% Manually plot the 5th bar filled
bar(5, betaW45Harms(5), 'FaceColor', betaW45Color, 'EdgeColor', betaW45Color);

% Get the number of bars
numBars = length(b45.YData);

% Loop through each bar and modify its appearance
for i = 1:numBars
    % Get the X and Y data of the current bar
    xVertices = get(b45, 'XData');
    % Calculate the X position for the current bar
    xPos = xVertices(i);
    
    % Draw top and side edges manually
    plot([xPos-0.4, xPos-0.4, xPos+0.4, xPos+0.4], ...
         [0, betaW45Harms(i), betaW45Harms(i), 0], ...
         'Color', betaW45Color, 'LineWidth', lineWidth);
end

hold off;                 
ylabel("Harm");
ylim([0,27]);
xticks([]);     
xticklabels([]);
if writeXTicks == 1
    xticks(1:5);
    xticklabels(runTypes);
end



figure;

% Create the bar plot
b376 = bar(betaW376Harms, 'LineWidth', lineWidth);
b376.FaceColor = 'none';  % Make all bars empty
b376.EdgeColor = 'none';  % Remove the default edges
hold on;

% Manually plot the 5th bar filled
bar(5, betaW376Harms(5), 'FaceColor', betaW376Color, 'EdgeColor', betaW376Color);

% Get the number of bars
numBars = length(b376.YData);

% Loop through each bar and modify its appearance
for i = 1:numBars
    % Get the X and Y data of the current bar
    xVertices = get(b376, 'XData');
    % Calculate the X position for the current bar
    xPos = xVertices(i);
    
    % Draw top and side edges manually
    plot([xPos-0.4, xPos-0.4, xPos+0.4, xPos+0.4], ...
         [0, betaW376Harms(i), betaW376Harms(i), 0], ...
         'Color', betaW376Color, 'LineWidth', lineWidth);
end
hold off;                 
ylabel("Harm");
ylim([0,27]);
xticks([]);     
xticklabels([]);   
if writeXTicks == 1
    xticks(1:5);
    xticklabels(runTypes);
end

figure;

% Create the bar plot
b3 = bar(betaW3Harms, 'LineWidth', lineWidth);
b3.FaceColor = 'none';  % Make all bars empty
b3.EdgeColor = 'none';  % Remove the default edges
hold on;

% Manually plot the 5th bar filled
bar(5, betaW3Harms(5), 'FaceColor', betaW3Color, 'EdgeColor', betaW3Color);

% Get the number of bars
numBars = length(b3.YData);

% Loop through each bar and modify its appearance
for i = 1:numBars
    % Get the X and Y data of the current bar
    xVertices = get(b3, 'XData');
    % Calculate the X position for the current bar
    xPos = xVertices(i);
    
    % Draw top and side edges manually
    plot([xPos-0.4, xPos-0.4, xPos+0.4, xPos+0.4], ...
         [0, betaW3Harms(i), betaW3Harms(i), 0], ...
         'Color', betaW3Color, 'LineWidth', lineWidth);
end

ylabel("Harm");
ylim([0,27]);
xticks([]);     
xticklabels([]);
if writeXTicks == 1
    xticks(1:5);
    xticklabels(runTypes);
end


combinedFig = figure('Position', [100, 100, 500, 125]);
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
    set(gca, 'FontSize', 12); % Changes font size for all text in the axes
end


fig = gcf; % Get current figure handle
fig.Position(3:4) = [1200 350]; % Set width and height in pixels


if getLegend == 1
    figure;

    % Create invisible plots with the desired styles
    plot(nan, nan, 'w'); 
    hold on;
    plot(nan, nan, 'w');
    plot(nan, nan, 'w');
    plot(nan, nan, 'w');
    plot(nan, nan, 'w');
    legendNames = {'No Intervention Policy', 'Optimal Fixed Policy', 'Optimal Stepwise Policy', 'Optimal Continuous Policy', 'Maximum Employment Reduction'};
    
    % Manually create the legend
    legend(legendNames);
    
    % Adjust the legend position or other properties if needed
    legend('Location', 'best');
end
