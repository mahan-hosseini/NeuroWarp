%% Plot comparison histogram plots of Marginals for N2pc and P3 separately
% Uses extenal nhist function
% ==> https://de.mathworks.com/matlabcentral/fileexchange/27388-plot-and-compare-histograms-pretty-by-default

function compareMarginals(filepath, separate)

% convert to char
filepath = char(filepath);

% load & check if plotpath is a dir for saving plots
load([filepath 'Marginal Distributions.mat'])
plotpath = [filepath 'Plots/Latency_correlations/'];
if ~exist(plotpath,'dir')
    mkdir(plotpath)
end

% If you want to plot histograms on the same plot as pdfs
if ~separate
    
    % N2pc
    fig = figure;
    colormap([0, 0, 0;1, 0, 0]) % to have one black and one red line
    nhist(N2pc_marginals, 'colors', 'colormap')
    fig.Children(2).XLabel.String = 'DTW Area';
    fig.Children(2).YLabel.String = 'Probability';
    fig.Children(2).Title.String = 'N2pc';
    fig.Children(2).FontSize = 14;
    saveas(gcf, [plotpath 'N2pc Marginal Comparison - Same Plot.jpg'])
    
    % P3
    fig = figure;
    colormap([0, 0, 0;1, 0, 0]) % to have one black and one red line
    nhist(P3_marginals, 'colors', 'colormap')
    fig.Children(2).XLabel.String = 'DTW Area';
    fig.Children(2).YLabel.String = 'Probability';
    fig.Children(2).Title.String = 'P3';
    fig.Children(2).FontSize = 14;
    saveas(gcf, [plotpath 'P3 Marginal Comparison - Same Plot.jpg'])
    

% If you want to plot histograms on separate plots with same number of bins    
elseif separate
    
    % N2pc
    fig = figure;
    nhist(N2pc_marginals,'color',rgb('lavender'),'separate','samebins','maxbins',50)
    fig.Children(1).YLabel = []; fig.Children(2).YLabel = [];
    fig.Children(2).Title.String = {'N2pc', 'Correct'};
    hY = suplabel('Frequency', 'y', [.12 .08 .84 .84]);
    hX = suplabel('DTW Area', 'x', [.08 .12 .84 .84]);
    hY.FontSize = 12; hX.FontSize = 12; 
    saveas(gcf, [plotpath 'N2pc Marginal Comparison - Separate Plots.jpg'])
    
    % P3
    fig = figure;
    nhist(P3_marginals,'color',rgb('lavender'),'separate','samebins','maxbins',50)
    fig.Children(1).YLabel = []; fig.Children(2).YLabel = [];
    fig.Children(2).Title.String = {'P3', 'Correct'};
    hY = suplabel('Frequency', 'y', [.12 .08 .84 .84]);
    hX = suplabel('DTW Area', 'x', [.08 .12 .84 .84]);
    hY.FontSize = 12; hX.FontSize = 12;
    saveas(gcf, [plotpath 'P3 Marginal Comparison - Separate Plots.jpg'])
end
end