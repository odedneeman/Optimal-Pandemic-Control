% This is Figure 2, the first plot.

% We have 2 panels, each with 3 plots. A planel for betaW and a panel for
% betaN.
% We plot death toll by employment reduction, GDP loss by emp reduction,
% and the harm by emp reduction, all for fixed policies- i.e. fixed emp
% reductions.




close all;

figureType = 1; 
fileName = "narrative_plot_range_betaN_logScaleDeathToll_13.6.png";
writeFile = 0;
getLegend = 0;
paramType = "betaW";
asterisk_legend = 0;
lineWidth = 2;
fontSize = 15;
xTicks = 0:10:50;
if figureType == 1

    load("ContinuousTime\statsPerN_Benchmark.mat");
    
    statsPerN = real(statsPerN); %technical matter
    nLevels = statsPerN(1, :);
    deathToll = statsPerN(2, :);
    GDPLoss = statsPerN(3, :);
    harm = statsPerN(4, :);
    span = statsPerN(5, :);
    
    nToDraw = {0.68, 0.8, 0.9, 1};
    marker_sizes = {0.2, 1/3, 2/3, 1};  % Sizes for empty, 1/3 full, 2/3 full, full
    
    nReduction = 100*(1 - nLevels); 

    figure(1)
    plot(nReduction, deathToll, 'k', 'LineWidth', lineWidth);
    hold on
    xlim([0, 32]);
    xlabel("Employment Reduction (in %)");
    ylabel({"Death Toll", "(per 100k)"});
    %set(gca, 'YScale', 'log');
    set(gca, 'FontSize', fontSize);
    %set(gca, 'FontName', 'Palatino Linotype', 'FontSize', 10);
    set(gca, 'XTick', xTicks); 
    if paramType == "betaN"
        yticks = 0:50:350;
        
        ylabels = {'0', '', '100', '', '200', '', '300', ''};
        
        set(gca, 'YTick', yticks, 'YTickLabel', ylabels);
    end

    figure(2)
    plot(nReduction, GDPLoss, 'k', 'LineWidth', lineWidth);
    hold on
    xlim([0, 32]);
    xlabel("Employment Reduction (in %)");
    ylabel({"GDP Loss", "(% of annual GDP)"});
    yticks = 0:5:50;
    ylabels = {'0', '', '10', '', '20', '', '30', '', '40', '', '50'};
    set(gca, 'YTick', yticks, 'YTickLabel', ylabels);
    ylim([0, 50]);
    set(gca, 'FontSize', fontSize);
    %set(gca, 'FontName', 'Palatino Linotype', 'FontSize', 10);
    set(gca, 'XTick', xTicks); 

    figure(3)
    plot(nReduction, harm, 'k', 'LineWidth', lineWidth);
    hold on
    set(gca, 'FontSize', fontSize);
    %set(gca, 'FontName', 'Palatino Linotype', 'FontSize', 10);
    set(gca, 'XTick', xTicks);

    load("singleSummaryStats\benchmarkStats.mat")

    x_point = benchmarkStats(1);
    y_point = benchmarkStats(2);
    benchmarkHarm = benchmarkStats(4);

    
    figure(3)
    legendNames2 = {"Benchmark"};
    legendNames3 = {};
    xlabel("Employment Reduction (in %)");
    ylabel({"Welfare loss (% of equivalent", "consumption loss)"});
    xlim([0, 32]);
    ylim([0, 30]);
    if paramType == "betaN"
        yticks = 0:2.5:30;
        
        ylabels = {'0', '', '5', '', '10', '', '15', '', '20', '', '25', '', '30'};
        
        %set(gca, 'YTick', yticks, 'YTickLabel', ylabels);
    end

    %betaWColors = {'b', [0.93, 0.51, 0.93], [1, 0.64, 0], 'r'};
    betaWColors = {[47, 151, 193]/255, 'r'};
    %betaNColors = {'r', [1, 0.64, 0], 'b'};
    betaNColors = {[229, 99, 153]/255, [209, 123, 15]/255};
    %betaWList = {0.3, 0.35, 0.4, 0.45};
    betaWList = {0.3, 0.45};
    betaNList = {0.6, 0.4};
    paramList = betaWList;
    colors = betaWColors;
    if paramType == "betaN"
        paramList = betaNList;
        colors = betaNColors;
    end
    


    for i = 1:numel(paramList)
        paramValue = paramList{i};
        load("ContinuousTime\statsPerN_" +  paramType + "_" + paramValue + ".mat");
        nLevels = statsPerParam(1, :);
        deathToll = statsPerParam(2, :);
        GDPLoss = statsPerParam(3, :);
        harm = statsPerParam(4, :);
        span = statsPerParam(5, :);
        
        nReduction = 100*(1 - nLevels);
        figure(1)
        hold on
        plot(nReduction, deathToll, 'Color', colors{i}, 'LineWidth', lineWidth);
        if paramType == "betaW"
            ylim([0, 450]);
        end
        if paramType == "betaN"
            ylim([0, 350]);
        end
        %legendNames1{end+1} = paramType + " = " + paramValue;
            
        figure(2)
        hold on
        plot(nReduction, GDPLoss, 'Color', colors{i}, 'LineWidth', lineWidth);
        
        figure(3)
        hold on
        plot(nReduction, harm, 'Color', colors{i}, 'LineWidth', lineWidth);
        legendNames2{end+1} = paramType + " = " + paramValue;


    end

    

    if getLegend==1
        figure(1)
        fig = gcf;
        legendNames1 = {"Low Disease Transmission", "Medium Disease Transmission", "High Disease Transmission"};
        legend(legendNames1, 'Location', 'southoutside', 'fontsize', 8);
    end

    figure(3)
    if getLegend==1
        legend(legendNames2);
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
    set(gca, 'FontSize', 12); % Changes font size for all text in the axes
end


fig = gcf; % Get current figure handle
fig.Position(3:4) = [1000 300]; % Set width and height in pixels
if writeFile == 1
    print(fileName, '-dpng', '-r600'); % 600 dpi resolution
end
end

figure;

% Create invisible plots with the desired styles
plot(NaN, NaN, 'Color', colors{1}, 'LineWidth', lineWidth); 
hold on;
plot(NaN, NaN, 'Color', 'k', 'LineWidth', lineWidth); 
plot(NaN, NaN, 'Color', colors{2}, 'LineWidth', lineWidth); 

legendNames = {'Low  Basic  Disease  Transimission', 'Medium  Basic  Disease  Transimission', 'High  Basic  Disease  Transmission'};

if paramType == "betaN"
    legendNames = {'High  Employment  Impact  Level', 'Medium  Employment  Impact  Level', 'Low  Employment  Impact  Level'};
end
% Manually create the legend
legend(legendNames);

% Adjust the legend position or other properties if needed
legend('Location', 'best');


