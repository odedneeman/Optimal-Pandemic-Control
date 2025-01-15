% Figure S2: time plots.
close all;

lineWidth = 2;
betaW3Color = [47, 151, 193]/255;
betaW376Color = 'k';
betaW45Color  = 'r';

colorNames = {betaW3Color, betaW376Color, betaW45Color};

betaWList = {0.3, 0.376, 0.45};
betaWRunNames = {};
eulerConst = 0.577216;

for i = 1:numel(betaWList)
    beta = betaWList{i};
    betaWRunNames{end+1} = "DiscreteTime\DiscreteSensitivityResultsWithEndo\beta_W_" + beta + ".mat";
end

runNames = betaWRunNames;
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
    param = betaWList{run};
    load(runNames{run});
    tV = round(muTv - eulerConst * sigmaTv);
    [betaListVac, xListVac, nListVac] = getVaccineState(xList, betaList, nList, tV, simulationDt, finalStep);
    figure(1);
    employmentReduction = 100*(1 - nListVac);
    plot(timeList, employmentReduction, 'Color', colorNames{run}, 'LineWidth', lineWidth); % Employment
    EffectiveR = (betaListVac / gamma) .* xListVac(1, :);
    EffectiveR(tV*100:end) = 0; % Forcing R-Effective to be 0 after vaccine arrival
    figure(2);
    plot(timeList(1:tV*100-1), EffectiveR(1:tV*100-1), 'Color', colorNames{run}, 'LineWidth', lineWidth); % R-Effective until the vaccine
    figure(3);
    plot(timeList, 1 - xListVac(1, :), 'Color', colorNames{run}, 'LineWidth', lineWidth); % 1-S
    figure(4);
    dDot = delta * theta * xListVac(4, :);
    plot(timeList, dDot * 100000, 'Color', colorNames{run}, 'LineWidth', lineWidth); % dDot per 100k
    figure(5);
    load("accumHarm\sensitivity\betaW_" + param + ".mat");
    plot(timeList, accumLoss, 'Color', colorNames{run}, 'LineWidth', lineWidth);
end
figure(1)
title("Employment Reduction");
xlim([0,finalTime]);
ylim([0,26]);
legend(runNames);

figure(2)
title("Effective R");
xlim([0,finalTime]);
ylim([0.67, 1.2]);
plot(timeList, ones(size(timeList)), 'k:', 'LineWidth', 1.5); % Adding a line for 1
legend([runNames, "R-Effective=1"]);

figure(3)
title("Accumulated Infections");
xlim([0,finalTime]);
ylim([0,0.36]);
legend(runNames);

figure(4)
title("Death Rate per 100K");
xlim([0,finalTime]);
ylim([10^-3, 1])
legend(runNames);
set(gca, 'YScale', 'log');
ytickformat('%,.0f');

figure(5)
title("Accumulated Welfare Loss");
xlim([0,finalTime]);
ylim([0, 24]);
legend(runNames);
ytickformat('%,.0f');

% Create a new figure
combinedFig = figure('Position', [0, 0, 500, 125]);

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
    set(ax, 'FontSize', 10);

    % Close the original figure to avoid cluttering the workspace
    close(figHandle);
end
%legend(legendNames);
fig = gcf; % Get current figure handle
fig.Position(3:4) = [1500 250]; % Set width and height in pixels
%subplot(1, 5, 4);
%yticks(10.^(0:3));