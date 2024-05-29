% Input Parameters
% ----------------
% 1) analysis_design
%    => Character array or string informing about your analysis design
%    => Has to be either "within" or "between", based on whether you have
%       a within- or between-subjects design
% 2) query
%    => Query time series (matrix w. dimensions: datapoints x subjects)
% 3) reference
%    => Reference time series (matrix w. dimensions: datapoints x subjects)
% 4) name_query
%    => Name of the query time series
% 5) name_reference
%    => Name of the reference time series
% 6) units
%    => In what units was your data measured? 
%       - For microvolts enter '\muV'
% 7) sampling_rate
%    => What was the sampling rate of your signal? In Hertz
%    => Used to compute the latency difference in milliseconds
% 8) savepath
%    => Path with files
% 9) permutations
%    => Number of permutations to adopt in statistical testing 
%       - We recommend at least 10000

% Note
% ----
% 1) This whole function is in line with the proposition of 
%    Zoumpoulaki et al. (2015) 
%    Latency as a region contrast: Measuring ERP latency differences with 
%    Dynamic Time Warping
% 2) Your full time series are analysed. 
%    If you would like to assess only a specific interval of your time 
%    series, use indexing before calling this function.

function neurowarp_latency_difference(analysis_design, query, reference, name_query, name_reference, units, sampling_rate, savepath, permutations)
%%       LATENCY DIFFERENCES WITH DYNAMIC TIME WARPING (IN MS)

%%                      SOME PREPARATION
% Check if vars were given as a string and not char, if not convert
if isstring(analysis_design)
    analysis_design = char(analysis_design);
end
if isstring(units)
    units = char(units);
end
if isstring(name_query)
    name_query = char(name_query);
end
if isstring(name_reference)
    name_reference = char(name_reference);
end
if isstring(savepath)
    savepath = char(savepath);
end
if (savepath(end) ~= '\') | (savepath(end) ~= '/')
    savepath = strcat(savepath, '/');
end

if ~ (strcmp(analysis_design, 'between') | strcmp(analysis_design, 'within'))
    disp('analysis_design must be either "between" or "within"! Fix & re-run')
    return
end
% Quick sanity check if series are identical for within and have same
% length (can have different number of subjects) if between
if strcmp(analysis_design, 'within')
    if  size(query) ~= size(reference)
        disp('Note - query & reference don''t have identical shapes! Fix & re-run!')
        return
    end
else % has to be 'between' now (see check above)
    if size(query, 1) ~= size(reference, 1) 
        disp('Note - query & reference do not have same datapoint-length! Fix & re-run!')
        return
    end
end

% Get averages
%   ==> plot_avgs are unstandardized for plotting
%   ==> z_avg immediately before running DTW
%   ==> QUERY is on X-AXIS 
%   ==> REFERENCE is on Y-AXIS

plot_avg_query = mean(query,2);  % unstandardized for plotting WP
plot_avg_reference = mean(reference,2);
z_avg_query = zscore(plot_avg_query); % standardized for analysis
z_avg_reference = zscore(plot_avg_reference); 

%       ---- MAIN DTW COMPUTATION HAPPENS HERE ----
% Outputs IX & IY == indices of x & y coordinates of warping path!
%   ORIGINAL FORM: [dist, ix, iy] = dtw(x,y)
%   ==> QUERY IS ON X-AXIS OF WP, SO FIRST INPUT
%   ==> REFERENCE IS ON Y-AXIS OF WP, SO SECOND INPUT!
%   ==> WP of this will be below main diagonal if QUERY is 
%       LATER than REFERENCE
%   ==> I.e. AREA will be positive in this case because area = areaDIAG -
%       areaWP (and hence areaDIAG will be larger than areaWP if WP is below
%       DIAG [and QUERY is indeed LATER than REFERENCE!]
[~, i_query_x, i_reference_y] = dtw(z_avg_query, z_avg_reference); 

%Plot original signals and aligned ones
%   Using following logic for plotting (from Matlab's DTW page)
% [d,i1,i2] = dtw(a1,a2); % do the DTW
% a1w = a1(i1);           % warped time-series if x-indices (of WP) of first series
% a2w = a2(i2);           % same logic for second series
warped_query = z_avg_query(i_query_x);
warped_reference = z_avg_reference(i_reference_y);

% plotting the first two plots (averages) 
figure;
set(gcf,'Position', [ 0        0        1280         907]);
subplot(4,1,1)
A_p(1) = plot(z_avg_query); 
ax(1) = gca;
subplot(4,1,2)
A_p(2) = plot(z_avg_reference);
ax(2) = gca;
% aligned series are longer, so we need a different x-axes for these
x_ax_aligned = 1:length(warped_query);
subplot(4,1,3)
A_p(3) = plot(x_ax_aligned, warped_query);
ax(3) = gca;
subplot(4,1,4)
A_p(4) = plot(x_ax_aligned, warped_reference);
ax(4) = gca;
hold on
% ylims
ylim_Real(1:2) = 0; ylim_Align(1:2) = 0;
ylim_Real(1) = min(-1, min(min([z_avg_query, z_avg_reference]))); 
ylim_Real(2) = max(1, max(max([z_avg_query, z_avg_reference])));
ylim_Align(1) = min(-1, min(min([warped_query, warped_reference])));
ylim_Align(2) =  max(1, max(max([warped_query, warped_reference])));
for i = 1:length(ax)
    ax(i).FontSize = 15;
    switch i 
        case 1
            A_p(i).Color = 'k';
            ax(i).YLim = ylim_Real;
        case 2
            A_p(i).Color = 'r';
            ax(i).XLabel.String = 'Time (ms)'; ax(i).XLabel.FontSize = 20;
            ax(i).YLim = ylim_Real;
        case 3
            A_p(i).Color = 'k';
            ax(i).YLim = ylim_Align;
        case 4
            A_p(i).Color = 'r';
            ax(i).XLabel.String = 'Datapoints'; ax(i).XLabel.FontSize = 20;
            ax(i).YLim = ylim_Align;
    end
A_p(i).LineWidth = 2;      
end
% suplabel & save
[~, h1] = suplabel('Alignment', 't', [.08 .12 .84 .84] );
h1.FontSize = 25;
[~, h2] = suplabel([units ' (z-scored)'], 'y', [.12 .08 .84 .84]);
h2.FontSize = 25;
saveas(gcf, [savepath 'Alignment.jpg'])
hold off

% Latency in MS is distribution of point-wise distances (distance = ix - iy) / (sampling rate (Hz) / 1000)
distance = i_query_x - i_reference_y;
latency = median(distance);
lat_in_ms = latency / (sampling_rate / 1000); 
% this is as proposed in Zoumpoulaki et al. (2015)

fprintf([' \n *** Latency Difference: ' num2str(lat_in_ms) 'ms! *** \n'])

% Plot warping paths compared to diagonal
plotWP(i_query_x, i_reference_y, plot_avg_query, plot_avg_reference, name_query, name_reference, units, lat_in_ms, savepath);

%%             PERMUTATION TEST ON AREA MEASURE (have this for median(distance), i.e. lat_in_ms measure, too)
%       GET OBSERVED AREA
DIAG = 1:max(i_query_x); % give me the diagonal
areaDIAG = trapz(DIAG); % area under diagonal
areaWP = trapz(i_query_x, i_reference_y);
area = (areaDIAG - areaWP) / areaDIAG; % area between diagonal and WP (whether this value is positive or negative informs about temporal order)

%       GET PERMUTED AREA
original_query = query;
original_reference = reference;

% Initialise
perm_x = cell(1,permutations); perm_y = cell(1, permutations); perm_area = cell(1, permutations);
datapoints = size(query, 1);
subjects = size(query, 2);
% Initialize permutation variables
for perms = 1:permutations
    % For a within-subject contrast, do the coin flipping for permuting
    % class labels
    if analysis_design == "within"
        perm_query(1:datapoints, 1:subjects) = 0;
        perm_reference(1:datapoints, 1:subjects) = 0;
        for s = 1:subjects
            coin = randi(2, subjects, 1);  % flip the coin
            if coin(s) == 1
                perm_query(:, s) = original_query(:, s); % assign based on coin
                perm_reference(:, s) = original_reference(:, s);  
            elseif coin(s) == 2
                perm_query(:, s) = original_reference(:, s);
                perm_reference(:, s) = original_query(:, s);
            end
        end
    % For a between-subject contrast, randomly assign subjects to groups,
    % preserving the original group sizes of query & reference
    else
        subjects_query = size(query, 2);
        subjects_reference = size(reference, 2);
        subjects_total = subjects_query + subjects_reference;
        this_perm = randperm(subjects_total);
        data_total = cat(2,query,reference);
        perm_query = data_total(:, this_perm(1:subjects_query));
        perm_reference = data_total(:, this_perm(subjects_query+1:subjects_total));
    end
    % Get averages
    perm_query_avg = mean(perm_query, 2);
    perm_reference_avg = mean(perm_reference, 2);
    % z-score permutation averages
    perm_query_z_avg = zscore(perm_query_avg);
    perm_reference_z_avg = zscore(perm_reference_avg);
    
    % DTW
    [~, perm_x{perms}, perm_y{perms}] = dtw(perm_query_z_avg, perm_reference_z_avg);
    
    % DTW area
    perm_DIAG = 1:max(perm_x{perms}); % diagonal
    perm_areaDIAG = trapz(perm_DIAG); % area under diagonal
    perm_areaWP = trapz(perm_x{perms}, perm_y{perms}); % area under WP
    perm_area{perms} = (perm_areaDIAG - perm_areaWP) / perm_areaDIAG; % standardized area inbetween permuted WP & diag
end

permdist_areas = cell2mat(perm_area);
permdist_areas = sort(permdist_areas); % sort the non-abs ones first (for indexing thresholds)
sort_abs_permdist_areas = sort(abs(permdist_areas));
thresh_abs = round(length(permdist_areas)*0.95);
thresh_abs_val = sort_abs_permdist_areas(thresh_abs); % find the 5% significance threshold using sorted absolutes
thresh1 = -thresh_abs_val; % find the first index that corresponds to -abs_thresh (bc we start left (negative))
thresh2 = thresh_abs_val; % find the first index that corresponds to +abs_thresh (starting at right (largest) now)

% area being positive or negative indicates latency-difference direction
%   we are doing two-tailed testing
p = sum(abs(permdist_areas) >= abs(area)) / permutations;

% get a string informing about p-value
clear plot_p
if p < 0.0001
    plot_p = 'p < 0.0001';
else
    plot_p = ['p = ' num2str(p)];
end

fprintf(['\n *** ' plot_p ' *** \n'])

%   Plot figure of observed area-size & its permutation-distribution
figure;
set(gcf,'Position', [ 0        0        1280         907]);
h = histogram(permdist_areas);
hold on
h.FaceColor = '#c7fdb5'; % pale green
h.EdgeColor = '#e4cbff'; % pale lilac  %h.EdgeColor = 'none';
t1_pos= vline(thresh1); t2_pos = vline(thresh2);
t1_pos.LineWidth = 3; t1_pos.LineStyle = '-'; t1_pos.Color = '#ff796c'; % salmon
t2_pos.LineWidth = 3; t2_pos.LineStyle = '-'; t2_pos.Color = '#ff796c'; % salmon
obs_line = vline(area);
obs_line.LineWidth = 3; obs_line.LineStyle = '-'; obs_line.Color = '#8ab8fe'; % carolina blue
ax = gca;
ax.FontSize = 17;
yticklabels_num = ax.YTick / permutations;
yticklabels_str = string(yticklabels_num);
ax.YTickLabels = yticklabels_str;
ax.YLabel.String = 'Proportion';
ax.XLabel.String = 'DTW - Area';
title({['Latency diff.: ' num2str(lat_in_ms) 'ms - ', plot_p]}, 'FontSize', 20)
saveas(gcf, [savepath 'Permutation Distribution.jpg'])
hold off
 

%   Plot same figure again but now plot the absolute permdiff values!
figure;
set(gcf,'Position', [ 0        0        1280         907]);
h = histogram(abs(permdist_areas));
hold on
h.FaceColor = '#c7fdb5'; % pale green
h.EdgeColor = '#e4cbff'; % pale lilac  %h.EdgeColor = 'none';
t1_pos= vline(thresh_abs_val); 
t1_pos.LineWidth = 3; t1_pos.LineStyle = '-'; t1_pos.Color = '#ff796c'; % salmon
obs_line = vline(abs(area));
obs_line.LineWidth = 3; obs_line.LineStyle = '-'; obs_line.Color = '#8ab8fe'; % carolina blue
ax = gca;
ax.FontSize = 17;
yticklabels_num = ax.YTick / permutations;
yticklabels_str = string(yticklabels_num);
ax.YTickLabels = yticklabels_str;
ax.YLabel.String = 'Proportion';
ax.XLabel.String = 'DTW - Area (absolute)';
title({['Absolute Areas - Latency diff.: ' num2str(lat_in_ms) 'ms - ', plot_p]}, 'FontSize', 20)
saveas(gcf, [savepath 'Permutation Distribution (absolutes).jpg'])
hold off

%   Save permutation-variables to file
save([savepath 'DTW Permutation Test - ' analysis_design ' subjects.mat'], 'area', 'p', 'thresh1', 'thresh2', 'perm_x', 'perm_y', 'permdist_areas')
clear area p thresh
end

%% PLOT WARPING PATH (with non-standardized GAs on Axes)!
function plotWP(i_query_x, i_reference_y, plot_avg_query, plot_avg_reference, name_query, name_reference, units, lat_in_ms, savepath)
%% Plot warping path compared to diagonal
figure;
set(gcf,'Position', [ 0        0        1280         907]);
% WP
subplot(4,4, [2:4 6:8 10:12])
wp_fontsize = 17;
% Set some axes stuff
ax1 = gca;
ax1.XLim(2) = i_query_x(end); ax1.YLim(2) = i_reference_y(end);
xmax = ax1.XLim(2);
hold on
% plot main diag
WPLineWidth = 3;
diag = plot(1:xmax, 1:xmax, 'r');
diag.Color = '#ff796c'; % salmon
diag.LineWidth = WPLineWidth;
% plot WP, 
p = plot(i_query_x, i_reference_y);
p.Color = '#8ab8fe'; % carolina blue
p.LineWidth = WPLineWidth;
% fill area between WP & main diagonal
f = fill(i_query_x, i_reference_y, rgb('pale green')); % pale green
f.FaceAlpha = 0.5;
f.EdgeAlpha = 0; % make this invisible
% this needs to be after plotting WP
ax1.XTickLabel = []; ax1.YTickLabel = []; % no tick labels for WP (but for GAs, below)

% Plot query on x-axis
subplot(4,4, 14:16)
pR = plot(plot_avg_query);
pR.Color = 'r'; pR.LineWidth = 3;
ylabel(units)
ax2 = gca;

% Plot reference on y-axis
subplot(4,4, [1 5 9])
pQ = plot(plot_avg_reference);
pQ.Color = 'k'; pQ.LineWidth = 3;
ylabel(units)
camroll(90)
ax3 = gca;
ax3.XAxisLocation = 'top';

ylimmax = max(max([plot_avg_query, plot_avg_reference])) + .5;
ax2.YLim(2) = ylimmax; ax3.YLim(2) = ax2.YLim(2);
ax2.XTick = ax1.XTick; ax3.XTick = ax2.XTick;
ax2.XLim(2) = length(plot_avg_query); ax3.XLim(2) = ax2.XLim(2);
ax2.FontSize = wp_fontsize; ax3.FontSize = wp_fontsize;
% add some (axis-title) text in plots 
t1 = text(ax2, ax2.XLim(2) - 30, ylimmax + 0.28, name_query);
t2 = text(ax3, ax3.XLim(2) + 15, ylimmax - 0.5, name_reference);
t1.FontSize = 20; t2.FontSize = t1.FontSize;
t1.FontAngle = 'italic'; t2.FontAngle = t1.FontAngle;
t1.FontWeight = 'bold'; t2.FontWeight = t1.FontWeight;
hold off

% set superlabels
hold on
[~, h1] = suplabel(['Latency difference: ' num2str(lat_in_ms) 'ms'], 't', [.2 .11 .85 .85]);
h1.FontSize = 22;
[~, h2] = suplabel('Time (ms)', 'x', [.08 .10 .85 .85]);
h2.FontSize = 20; 
[~, h3] = suplabel('Time (ms)', 'y');
h3.FontSize = h2.FontSize;
hold off

% Save it
saveas(gcf, [savepath 'WarpPath.jpg'])
 
end