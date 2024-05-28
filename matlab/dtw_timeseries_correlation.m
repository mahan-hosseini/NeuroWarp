% Conceptual Approach 
% -------------------
% To assess latency correlations between 2 time series:
%   1) Get Single Subject (or group) time series
%   2) Get N bootstrap time series (w. same bootstrapping 
%       subject-distributions for the two time series)
%   3) DTW of 2) with 1) using AREA measure!
%   4) Gives you two numbers for each bootstrap sample -> DTW result of 
%       time series 1 & 2 latency (lag)
%   5) Correlate these with Pearson and Spearman (latter preferred for
%       presentations)

% Input Parameters
% ----------------
% 1) series_1, series_2
%    => Timeseries as two matrices of shape: 
%       dataindices (e.g., timepoints) x subjects (or repetitions)
% 2) name_1, name_2
%    => Names of your timeseries
% 3) savepath
%    => Where to save output to
% 4) num_boots
%    => Number of bootstraps (we recommend at least 10000)
% 5) outlier
%    => 1: to exclude +-5SD outliers from your scatterplot
%    => 0: don't exclude outliers
% 6) try_to_fix_ylims
%    => 1: try to fix ylimits of marginals and include lines for easier
%          comparison 
%          - based on num_boots
%          - if you want to modify this yourself (e.g., hardcode lines &
%            maxima, change the if statement blocks for the two marginals)
%    => 0: let matlab handle ylimits of marginals (likely won't allow 
%          comparisons

% Note
% ----
% 1) Number of subjects has to be equal in both time series!
% 2) Your full time series are analysed. 
%    If you would like to assess only a specific interval of your time 
%    series, use indexing before calling this function.
% 3) We fix marginals' x-axis ranges to +-8 SD. This should be okay for 
%    most data-sets since we zscore DTW area distributions. However, if you
%    do not want to standardise like this (which will make your marginals 
%    harder to compare) or you want to implement different values, change 
%    the lines that use xlimits & ylimits (see line 180)

function dtw_timeseries_correlation(series_1, series_2, name_1, name_2, savepath, num_boots, outlier, try_to_fix_ylims)
%%                          PREPARATION

% first check if user adhered to equal sample size across time series
if size(series_1, 2) ~= size(series_2, 2)
    disp("Your time series have a different number of subjects - cancelling!")
    return
end

% then check if names were strings, have to be chars
if isstring(name_1)
    name_1 = char(name_1);
end
if isstring(name_2)
    name_2 = char(name_2);
end

% prepare path for plots
if savepath(end) ~= '/'
    savepath = strcat(savepath, '/');
end
plotpath = [savepath 'Plots/'];
if ~exist(plotpath,'dir')
    mkdir(plotpath)
end
if ~exist([savepath 'Variables'], 'dir')
    mkdir([savepath 'Variables'])
end

savediaryto = [savepath 'CorrelateAnalysis_Info.txt'];
diary(savediaryto)

% for saving marginals & stats to file
Marginals = struct(name_1, [], name_2, []);
Stats = struct();

%%                  RUN THE ANALYSIS 
% run DTW & bootstrap
[series_1_dtw_area, series_2_dtw_area] = run_bootstrap_dtw(series_1, series_2, num_boots);

%% OUTLIER REJECTION? 
if outlier
    % exclude outliers if +-5 SD (change in std_thresh line if wanted)

    series_1_std = std(series_1_dtw_area);
    series_2_std = std(series_2_dtw_area);
    std_thresh = 5;                 % threshold for outlier rejection!

    % *****************************************************************
    %                           IMPORTANT 
    % This might seem like a strange way of removing outliers 
    % ==> However it ensures that we check all 4 directions of the
    %   scatterplot for outliers, save indices as we do so, then just
    %   take those unique indices (removing duplicates) and (IMPORTANTLY)
    %   IN THE END removing these indices FROM BOTH DTW distributions
    % ==> THIS IS IMPORTANT BECAUSE DTW DISTRIBUTIONS HAVE TO STAY
    %   LINKED DUE TO THE BOOTSTRAP
    % *****************************************************************
    
    outlier_idx = []; 
    series_1_pos_outlier_idx = find(series_1_dtw_area > (mean(series_1_dtw_area) + (std_thresh*series_1_std)));
    if ~isempty(series_1_pos_outlier_idx)
        outlier_idx = vertcat(outlier_idx, series_1_pos_outlier_idx);
    end
    series_1_neg_outlier_idx = find(series_1_dtw_area < (mean(series_1_dtw_area) - (std_thresh*series_1_std)));
    if ~isempty(series_1_neg_outlier_idx)
        outlier_idx = vertcat(outlier_idx, series_1_neg_outlier_idx);
    end
    series_2_pos_outlier_idx = find(series_2_dtw_area > (mean(series_2_dtw_area) + (std_thresh*series_2_std)));
    if ~isempty(series_2_pos_outlier_idx)
        outlier_idx = vertcat(outlier_idx, series_2_pos_outlier_idx);
    end
    series_2_neg_outlier_idx = find(series_2_dtw_area < (mean(series_2_dtw_area) - (std_thresh*series_2_std)));
    if ~isempty(series_2_neg_outlier_idx)
        outlier_idx = vertcat(outlier_idx, series_2_neg_outlier_idx);
    end
    outlier_idx = unique(outlier_idx); % get rid of duplicates
    
    series_1_dtw_area(outlier_idx) = []; % remove outlier indices
    series_2_dtw_area(outlier_idx) = [];
end

%%         Save Var, Skew & Kurtosis Values of DTW Area Distributions
Stats = struct('series_1_var', var(series_1_dtw_area), 'series_1_skew', skewness(series_1_dtw_area), 'series_1_kurt', kurtosis(series_1_dtw_area), ...
    'series_2_var', var(series_2_dtw_area), 'series_2_skew', skewness(series_2_dtw_area), 'series_2_kurt', kurtosis(series_2_dtw_area));
if outlier
    save([savepath 'Variables/' num2str(std_thresh) 'SD OutlierRejected Marginal Stats.mat'], 'Stats') % save marginal-stats to file
else
    save([savepath 'Variables/Marginal Stats.mat'], 'Stats')
end

%%              Standardise DTW Area Distributions 
% ==> So line of linear fit in scatter & marginal dists reflect
%   correlation values better (since correlations have standardisation
%   internally)
x = series_1_dtw_area;
y = series_2_dtw_area;

% standardize x & y
x = zscore(x);
y = zscore(y);   

z_x = x; % for saving standardised and original marginals later
z_y = y;

%% *************************************************************************
%                           correlation stuff
[pears_r, pears_p] = corr(x, y, 'Type', 'Pearson');
[spear_r, spear_p] = corr(x, y, 'Type', 'Spearman');
pears_str_r = num2str(round(pears_r, 2));
spear_str_r = num2str(round(spear_r, 2));
str_df = num2str(num_boots-2);
if pears_p < 0.0001
    pears_str_p = 'p < .0001';
else
    pears_p = round(pears_p,3);
    pears_str_p = ['p = ' num2str(pears_p)];
end
if spear_p < 0.0001
    spear_str_p = 'p < .0001';
else
    spear_p = round(spear_p,3);
    spear_str_p = ['p = ' num2str(spear_p)];
end

%% Plot everything in one big figure
figure;
set(gcf,'Position', [ 0        0        1280         907]);
% some common configs
scatter_color = rgb('lavender');
line_color = rgb('carolina blue');
ax_fontsize = 20;  
xlimits = [-8 8]; ylimits = [-8 8]; % when zscoring marginals
if try_to_fix_ylims
    ylimmax = num_boots / 8;
    ref_val_1 = round(ylimmax*0.25);
    ref_val_2 = round(ylimmax*0.5);
    ref_val_3 = round(ylimmax*0.75);
end
%% first, scatterplot
subplot(4,4,[2:4 6:8 10:12])
scatter(x,y, [], scatter_color);
hold on
ax = gca;
ax.FontSize = ax_fontsize;
ax.YLim = ylimits; % when zscoring marginals
ax.XLim = xlimits;
series_1_xlimits = ax.XLim;   % use the scatterplot's xlimits for series_1's marginal xlimits
series_2_xlimits = ax.YLim; % use the scatterplot's Y (!) limits for series_2's marginal X (!) limits
ax.XTickLabels = [];
ax.YTickLabels = [];
title({['DTW area-diff correlations of ' name_1 ' & ' name_2 '.'], ...
    ['Pearson: r(' str_df ') = ' pears_str_r ': ' pears_str_p ' & Spearman: r(' str_df ') = ' ...
    spear_str_r ': ' spear_str_p]}, 'fontsize', 20)
l = lsline;
l.LineWidth = 3;
l.Color = line_color;
%% second, series_2 marginal distribution
% initialise axes
subplot(4,4,[1 5 9])
ax_series_2 = gca;
ax_series_2.FontSize = ax_fontsize;
ax_series_2.XLim = series_2_xlimits;
camroll(90)
hold on
% reference lines on marginal & fix yticks
if try_to_fix_ylims
    ax_series_2.YLim(2) = ylimmax;
    reflines = hline([ref_val_1 ref_val_2 ref_val_3], {'k', 'k', 'k'}); % for comparison
    for r = 1:length(reflines)
        reflines(r).LineWidth = 0.2;
        reflines(r).Color = rgb('teal');
        reflines(r).LineStyle = '--';
    end
    ax_series_2.YTick = [ref_val_1 ref_val_2 ref_val_3];
    ax_series_2.YTickLabels = {num2str(ref_val_1), num2str(ref_val_2), num2str(ref_val_3)};
    ax_series_2.YTickLabelRotation = 45;
end
% histogram
h_series_2 = histogram(y, 100);
h_series_2.FaceColor = scatter_color;
h_series_2.EdgeColor = line_color;
set(ax_series_2,'xaxislocation', 'top') % put xticklabels on top
% Instead of including var/skew/kurt of marginals in the title I will save
% them to a file
%     % variance, skewness & kurtosis                                       
%     var_series_2 = var(y);
%     var_series_2 = round(var_series_2,3);
%     skew_series_2 = skewness(y); % compute skewness
%     skew_series_2 = round(skew_series_2,3);
%     kurt_series_2 = kurtosis(y); % compute kurtosis
%     kurt_series_2 = kurt_series_2 - 3; % standardize it so if kurt == normal, it's 0
%     kurt_series_2 = round(kurt_series_2,3);
%     ax_series_2.XLabel.String = {'series_2 Marginal Distribution', ['Var: ' num2str(var_series_2) ' | Skew: ' num2str(skew_series_2) ' | Kurt: ' num2str(kurt_series_2)]};
ax_series_2.XLabel.String = [name_2 ' Marginal Distribution'];
ax_series_2.XLabel.FontSize = 17; % "Title" as XLabel because of axis-rotation
ax_series_2.XLabel.FontWeight = 'bold';
hold off
%% third, series_1 marginal distribution
subplot(4,4,14:16)
% initialise axes
ax_series_1 = gca;
ax_series_1.FontSize = ax_fontsize;
ax_series_1.XLim = series_1_xlimits;
hold on
% reference lines & yticks
if try_to_fix_ylims
    ax_series_1.YLim(2) = ylimmax;
    reflines = hline([ref_val_1 ref_val_2 ref_val_3], {'k', 'k', 'k'}); 
    for r = 1:length(reflines)
        reflines(r).LineWidth = 0.2;
        reflines(r).Color = rgb('teal');
        reflines(r).LineStyle = '--';
    end
    ax_series_1.YTick = [ref_val_1 ref_val_2 ref_val_3];
    ax_series_1.YTickLabels = {num2str(ref_val_1), num2str(ref_val_2), num2str(ref_val_3)};
    ax_series_1.YAxisLocation = 'right';
end
% histogram
h_series_1 = histogram(x, 100);
h_series_1.FaceColor = scatter_color;
h_series_1.EdgeColor = line_color;
% variance, skewness & kurtosis
%     var_series_1 = var(x);       % compute variance
%     var_series_1 = round(var_series_1,3);
%     skew_series_1 = skewness(x); % compute skewness
%     skew_series_1 = round(skew_series_1,3);
%     kurt_series_1 = kurtosis(x); % compute kurtosis
%     kurt_series_1 = kurt_series_1 - 3; % standardize it so if kurt == normal, it's 0
%     kurt_series_1 = round(kurt_series_1,3);
%     title({'series_1 Marginal Distribution', ['Var: ' num2str(var_series_1) ' | Skew: ' num2str(skew_series_1) ' | Kurt: ' num2str(kurt_series_1)]}, 'FontSize', 15)
title([name_1 ' Marginal Distribution'], 'FontSize', 17)
hold off
%% Reposition axes to remove space between subplots, add labels & save
scatterPos = ax.Position;
series_2_Pos = ax_series_2.Position;
series_1_Pos = ax_series_1.Position;
series_2_Pos(1) = scatterPos(1) - series_2_Pos(3);
series_1_Pos(2) = scatterPos(2) - series_1_Pos(4);
ax_series_2.Position = series_2_Pos;
ax_series_1.Position = series_1_Pos;
% save it
if ~exist(plotpath, 'dir')
    mkdir(plotpath)
end
if outlier
    saveas(gcf,[plotpath 'DTW Latency Correlation - ' name_1 ' & ' name_2 ' - ' num2str(std_thresh) 'std outlier removed.jpg'])
else
    saveas(gcf,[plotpath 'DTW Latency Correlation - ' name_1 ' & ' name_2 '.jpg'])
end
hold off
% save marginals to struct
series_2_marginals = y;
series_1_marginals = x;
z_series_2_marginals = z_y;
z_series_1_marginals = z_x;

if outlier
    save([savepath 'Variables/Marginal Distributions - ' num2str(std_thresh) 'std outlier removed.mat'], 'series_2_marginals', 'series_1_marginals', 'z_series_2_marginals', 'z_series_1_marginals') % save marginals to file
else
    save([savepath 'Variables/Marginal Distributions.mat'], 'series_2_marginals', 'series_1_marginals', 'z_series_2_marginals', 'z_series_1_marginals')
end
diary off
end

%%                      RUN THE DTW BOOTSTRAP ANALYSIS
function [series_1_dtw_area, series_2_dtw_area] = ...
    run_bootstrap_dtw(series_1, series_2, num_boots)

% compute standardised (subjects/reps) averages of both time series
avg_series_1 = mean(series_1,2);
avg_series_2 = mean(series_2,2);
z_avg_series_1 = zscore(avg_series_1);
z_avg_series_2 = zscore(avg_series_2);

% bootstrap manually using randi & run dtw between bootstrap & observed
% z_avg
num_subs = size(series_1, 2);
num_bootstraps = num_boots;
series_1_dtw_area(num_bootstraps,1) = 0; series_2_dtw_area(num_bootstraps,1) = 0; % initialize stuff
% and loop
for n = 1:num_bootstraps
    
    disp(['boots num #' num2str(n)]) % info
    
    % get num_subj bootstrap samples (using randi == with duplicates)
    this_trap = randi(num_subs,num_subs,1); % max-val = num_subs, size = num_subs x 1
    boot_series_1 = series_1(:, this_trap);
    boot_series_2 = series_2(:, this_trap);
    
    avg_boot_series_1 = mean(boot_series_1,2);
    avg_boot_series_2 = mean(boot_series_2,2);
    
    z_avg_boot_series_1 = zscore(avg_boot_series_1);
    z_avg_boot_series_2 = zscore(avg_boot_series_2);
    
    % DTW between GA_boot & GA_obs
    %   initialize
    ix = []; iy = []; distance = [];
    
    %       series_1
    %   IX & IY ARE indices of x & y coordinates of warping path!
    [~, ix, iy] = dtw(z_avg_boot_series_1, z_avg_series_1);
    DIAG = 1:max(ix); % get the diagonal
    areaDIAG = trapz(DIAG); % area under diag
    areaWP = trapz(ix, iy); % area under WP
    series_1_dtw_area(n) = (areaDIAG - areaWP) / areaDIAG; % to standardize, as in Zoumpoulaki et al 2015
    
    %   reset
    clear ix iy areaWP areaDIAG
    
    %       series_2
    [~, ix, iy] = dtw(z_avg_boot_series_2, z_avg_series_2);
    DIAG = 1:max(ix); % get the diagonal
    areaDIAG = trapz(DIAG); % area under diag
    areaWP = trapz(ix, iy); % area under WP
    series_2_dtw_area(n) = (areaDIAG - areaWP) / areaDIAG;
    
    %   If you are interested to investigate your data, comment out the 
    %   bit below and either write a save statement or add the two vars
    %   below to this local function's output (they are the bootstraped 
    %   averages - useful if you have funny shapes in DTW
    %   distributions you would like to understand better)
    % all_boots_series_1s{n} = z_avg_boot_series_1;
    % all_boots_series_2s{n} = z_avg_boot_series_2;
end
end