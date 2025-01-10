% Figure S1: The endogenous response

close all;


lineWidth = 1.5;
endoColor = 'b';
policyColor = 'k';
deathColor = 'r';
getLegend = 1;

equal_time1 = 73.29;
equal_time2 = 240.65;

equal_index1 = 100*equal_time1 + 1;
equal_index2 = 100*equal_time2 + 1;

load("ContinuousTime\ContinuousSensitivityResultsWithEndo\betaW_0.45_No_Intervention.mat");
nReduction = 100*(1 - nList);
figure(1);
plot(timeList, nReduction, 'Color', endoColor, 'LineWidth', lineWidth);
xlim([0, 635]);
hold on
plot(timeList, zeros(1,numel(nReduction)), 'Color', policyColor, 'LineWidth', lineWidth, 'LineStyle', '--');
xlabel('Time (days)');
ylabel('Employment Reduction (%)');
ylim([0, 27]);
%title('No Employment Policy');

deathRate1 = zeros(1, finalStep);
for curStep = 1 : finalStep
    curX = xList(:, curStep);
    curResolving = curX(4);
    dDot = delta * theta * curResolving;
    deathRate1(curStep) = dDot * 100000;
end


figure(2)
plot(timeList, deathRate1, 'Color', deathColor, 'LineWidth', lineWidth)
xlim([0, 635]);
ylim([0, 1.6]);
xlabel('Time (days)');
ylabel('Death Rate (per 100K)');
set(gca, 'YTick', 0:0.2:1.6);

load("DiscreteTime\DiscreteSensitivityResultsWithEndo\beta_N_0.53.mat");
figure(3);
nReduction = 100*(1 - nList);
plot(timeList(1:equal_index1), nReduction(1:equal_index1), 'Color', policyColor, 'LineWidth', lineWidth);
hold on
plot(timeList(equal_index1:equal_index2), nReduction(equal_index1:equal_index2), 'Color', policyColor, 'LineWidth', lineWidth, 'LineStyle', '--');
plot(timeList(equal_index2:end), nReduction(equal_index2:end), 'Color', policyColor, 'LineWidth', lineWidth);
xlim([0, 635]);
xlabel('Time (days)');
ylabel('Employment Reduction (%)');
ylim([0, 27]);
%title('Insufficiently restrictive policy');


load("DiscreteTime\DiscreteRobustnessResults\beta_W_benchmark_0.376_reality_0.45.mat");

finalStep = round(finalTime / simulationDt) + 1;
endoResponse = zeros(1, finalStep);
deathRate2 = zeros(1, finalStep);
for curStep = 1 : finalStep
    curTime = (curStep - 1) * simulationDt;
    curX = xList(:, curStep);
    curResolving = curX(4);
    curKappa = getTVKappa(curTime);
    dDot = delta * theta * curResolving;
    curEndo = curKappa * dDot;
    endoResponse(curStep) = 100*curEndo;
    deathRate2(curStep) = dDot * 100000;
end


figure(3);  
plot(timeList(1:equal_index1), endoResponse(1:equal_index1), 'Color', endoColor, 'LineWidth', lineWidth, 'LineStyle', '--');
plot(timeList(equal_index1:equal_index2), endoResponse(equal_index1:equal_index2), 'Color', endoColor, 'LineWidth', lineWidth);
plot(timeList(equal_index2:end), endoResponse(equal_index2:end), 'Color', endoColor, 'LineWidth', lineWidth, 'LineStyle', '--');


figure(4)
plot(timeList, deathRate2, 'Color', deathColor, 'LineWidth', lineWidth)
xlim([0, 635]);
xlabel('Time (days)');
ylabel('Death Rate (per 100K)');
ylim([0, 1.6]);
set(gca, 'YTick', 0:0.2:1.6);

combinedFig = figure('Position', [100, 100, 500, 125]);
rows = 2;
cols = 2;


% Loop through each subplot
for i = 1:4
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
fig.Position(3:4) = [800 800]; % Set width and height in pixels

if getLegend == 1
    % Create a sample plot
    figure;
    hold on;
    p1 = plot(nan, nan, 'Color', policyColor, 'LineWidth', lineWidth);
    p2 = plot(nan, nan, 'Color', endoColor, 'LineWidth', lineWidth);
    p3 = plot(nan, nan, 'w'); % Invisible plot
    %p4 = plot(nan, nan, 'Color', deathColor, 'LineWidth', lineWidth);
    
    
    % Create the custom legend
    legend([p1, p2, p3], {'Employment Reduction Policy', 'Endogenous Policy Reduction', 'Dashed - Latent Employment Reduction'});
    
    hold off;
end