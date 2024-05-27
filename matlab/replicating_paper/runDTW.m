function runDTW(experiment, filepath, dtwpath, permutations, smooth)
%%       DYNAMIC TIME WARPING
%%INPUTS
%   P3s             -- struct containing P3s at Fz, Cz, Pz
%   N2pcs           -- struct containing N2pcs
%   permutations    -- NUMBER OF PERMUTATIONS OF PERMUTATION TEST

%%                      NOTE
% This whole function is in line with the proposition of 
%   Zoumpoulaki et al. (2015) 
%   Latency as a region contrast: Measuring ERP latency differences with 
%   Dynamic Time Warping

%%                      SOME PREPARATION
%   --- NOTE FOR MYSELF: COR IS REF (Y-AXIS) & INT IS QUERY (X-AXIS!) ---
clc
allERPs = {'P3', 'N2pc'};
allP3chans = {'Fz', 'Cz', 'Pz'};
% Print some info
fprintf(['\n *** USING DYNAMIC TIME-WARPING & TWO-TAILED PERMUTATION TEST RUNNING ' num2str(permutations) ' PERMUTATIONS *** \n'])
fprintf('\n *** ADOPTING CONDITIONS: REFERENCE = CORRECT (y-axis) & QUERY = INTRUSIONS (x-axis) *** \n')
fprintf('\n *** HENCE IF DTW MEASURE IS POSITIVE, MEANS QUERY LATER THAN REFERENCE *** \n')

% convert strings to char
experiment = char(experiment);
filepath = char(filepath);
dtwpath = char(dtwpath);

% some adjustments if we are analyzing experiment 2
if strcmp(experiment, "Experiment 2")
    exp_string = '_Exp2';
    dtwpath = [dtwpath 'Experiment 2\'];
else
    exp_string = '';
end

% check if dtwpath exists
if ~isfolder(dtwpath)
    mkdir(dtwpath)
end

% load
if smooth
    loadstruct = load([filepath 'LPF_ERP_structs' exp_string]);
    P3s = loadstruct.LPF_P3s;
    N2pcs = loadstruct.LPF_N2pcs;
else
    loadstruct = load([filepath 'ERP_structs' exp_string]);
    P3s = loadstruct.P3s;
    N2pcs = loadstruct.N2pcs;
end

% define subjects after loading (because experiments have different subnums)
subjects = length(P3s);

for e = 1:length(allERPs)
% for e = 3 % EXPLORE CZ!
    ERP2do = allERPs{e};    % loop over P3 & N2pc analyses
    if strcmp(ERP2do, 'P3')
        for P3chan = 1:length(allP3chans) % loop over 3 P3 channels
            % print what we test
            fprintf(['\n ********* TESTING ' allP3chans{P3chan} ' ********* \n'])
            % Get full_SSERP _Cor & _Int
            clear analysis tmp1 tmp2 full_SSERP_Cor full_SSERP_Int timewindow dataindices savename
            tmp1 = cat(3,P3s.Correct); tmp2 = squeeze(tmp1(:,P3chan,:));
            full_SSERP_Cor = tmp2;
            clear tmp1 tmp2
            tmp1 = cat(3,P3s.Intrusion); tmp2 = squeeze(tmp1(:,P3chan,:));
            full_SSERP_Int = tmp2;
            % some more vars
            analysis = ['P3 at ' allP3chans{P3chan}];
            if ismember(P3chan, [1 2])
                timewindow = [302 800];
            else
                timewindow = [252 800];
            end
            % run it
            computeDTW(full_SSERP_Cor, full_SSERP_Int, timewindow, subjects, permutations, analysis, dtwpath) 
        end
    elseif strcmp(ERP2do, 'N2pc')
        fprintf('\n ********* TESTING N2pc ********* \n')
        % Get full_SSERP _Cor & _Int
        clear analysis  tmp1 tmp2 full_SSERP_Cor full_SSERP_Int timewindow dataindices savename
        tmp1 = cat(3,N2pcs.Correct); tmp2 = squeeze(tmp1);
        full_SSERP_Cor = tmp2;
        clear tmp1 tmp2
        tmp1 = cat(3,N2pcs.Intrusion); tmp2 = squeeze(tmp1);
        full_SSERP_Int = tmp2;
        % some more vars
        analysis = 'N2pc';
        timewindow = [152 400];
        % run it
        computeDTW(full_SSERP_Cor, full_SSERP_Int, timewindow, subjects, permutations, analysis, dtwpath)
    end
end
end

%%                      DTW ANALYSIS
function computeDTW(full_SSERP_Cor, full_SSERP_Int, timewindow, subjects, permutations, analysis, dtwpath)
% Get dataindices from timewindow
time = -98:2:800;
dataindices = find(time>=timewindow(1) & time<=timewindow(2));

% Extract time-window
for i = 1:subjects
    SSERP_Cor(:, i) = full_SSERP_Cor(dataindices, i);
    SSERP_Int(:, i) = full_SSERP_Int(dataindices, i);
end

% Get Grand Averages 
%   ==> plot_GAs are unstandardized for plotting
%   ==> zscore GAs immediately before running DTW
%   ==> QUERY INTRUSION (X-AXIS) is what we hypothesize to be later (INTRUSION)
%   ==> REFERENCE CORRECT (Y-AXIS) is what we hypothesize to be earlier (CORRECT)

plot_GA_RefCorY = mean(SSERP_Cor,2); % unstandardized for plotting WP
plot_GA_QuerIntX = mean(SSERP_Int,2);
GA_RefCorY = zscore(plot_GA_RefCorY); % standardized for analysis
GA_QuerIntX = zscore(plot_GA_QuerIntX);

%       ---- MAIN DTW COMPUTATION HAPPENS HERE ----
% Outputs IX & IY == indices of x & y coordinates of warping path!
%   ORIGINAL FORM: [dist, ix, iy] = dtw(x,y)
%   ==> QUERY/INTRUSION IS ON X-AXIS OF WP, SO FIRST INPUT
%   ==> REF/CORRECT IS ON Y-AXIS OF WP, SO SECOND INPUT!
%   ==> WP of this will be below main diagonal if intrusions (QUERY) are 
%       LATER than corrects (REFERENCE) 
%   ==> I.e. AREA will be positive in this case because area = areaDIAG -
%       areaWP (and hence areaDIAG will be larger than areaWP if WP is below
%       DIAG [and intrusions (Q) are indeed later than corrects(R)!]
[~, i_QuerIntX, i_RefCorY] = dtw(GA_QuerIntX, GA_RefCorY); 

%Plot original signals and aligned ones
%   Using following logic for plotting (from Matlab's DTW page)
% [d,i1,i2] = dtw(a1,a2); % do the DTW
% a1w = a1(i1);           % warped time-series if x-indices (of WP) of first ERPs time-series
% a2w = a2(i2);           % same logic for second ERP
GA_RefCorY_w = GA_RefCorY(i_RefCorY);
GA_QuerIntX_w = GA_QuerIntX(i_QuerIntX);


% dataindices are indices of this current ERP component, so we use it for
% plotting the first two plots (ERPs)
figure;
set(gcf,'Position', [ 0        0        1280         907]);
subplot(4,1,1)
A_p(1) = plot(dataindices, GA_RefCorY); 
ax(1) = gca;
subplot(4,1,2)
A_p(2) = plot(dataindices, GA_QuerIntX);
ax(2) = gca;
% aligned ERPs are longer, so we need a different x-axes for these
x_ax_aligned = 1:length(GA_RefCorY_w);
subplot(4,1,3)
A_p(3) = plot(x_ax_aligned, GA_RefCorY_w);
ax(3) = gca;
subplot(4,1,4)
A_p(4) = plot(x_ax_aligned, GA_QuerIntX_w);
ax(4) = gca;
hold on
% ylims
ylim_Real(1:2) = 0; ylim_Align(1:2) = 0;
ylim_Real(1) = min(-1, min(min([GA_RefCorY, GA_QuerIntX]))); 
ylim_Real(2) = max(1, max(max([GA_RefCorY, GA_QuerIntX])));
ylim_Align(1) = min(-1, min(min([GA_RefCorY_w, GA_QuerIntX_w])));
ylim_Align(2) =  max(1, max(max([GA_RefCorY_w, GA_QuerIntX_w])));
% xticklabels for realERPs
xticklab = zeros(1,length(ax(1).XTick));
for i = 1:length(xticklab)
    xticklab(i) = time(ax(1).XTick(i));
end
% axes & line stuff
for i = 1:length(ax)
    ax(i).FontSize = 15;
    switch i 
        case 1
            A_p(i).Color = 'k';
            ax(i).YLim = ylim_Real;
            ax(i).XTickLabels = xticklab;
            ax(i).XLim = [dataindices(1) dataindices(end)];
        case 2
            A_p(i).Color = 'r';
            ax(i).XLabel.String = 'Time (ms)'; ax(i).XLabel.FontSize = 20;
            ax(i).YLim = ylim_Real;
            ax(i).XTickLabels = xticklab;
            ax(i).XLim = [dataindices(1) dataindices(end)];
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
[~, h1] = suplabel(analysis, 't', [.08 .12 .84 .84] );
h1.FontSize = 25;
[~, h2] = suplabel('\muV (z-scored)', 'y', [.12 .08 .84 .84]);
h2.FontSize = 25;
if ~isfolder(dtwpath)
    mkdir(dtwpath)
end
saveas(gcf, [dtwpath analysis ' -- Alignment.jpg'])
hold off
close all




%   ==> QUERY INTRUSION (X-AXIS) is what we hypothesize to be later (INTRUSION)
%   ==> REFERENCE CORRECT (Y-AXIS) is what we hypothesize to be earlier (CORRECT)

% LATENCY BEING POSITIVE FITS, BECAUSE WE TAKE HORIZONTAL DISTANCE HERE & WP
% IS BELOW MAINDIAG BC INT IS L8R THAN COR, SO THE "TYPICAL" (i.e.,
% median) thing would be that ix is larger than iy (because our vertical
% line is going to the "right" from the main diagonal because ix is larger
% than iy because int was later)

% Latency in MS is distribution of point-wise distances (distance = ix - iy) / (sampling rate (Hz) / 1000)
distance = i_QuerIntX - i_RefCorY;
latency = median(distance);
lat_in_ms = latency / .5; % 500Hz
% this is as proposed in Zoumpoulaki et al. (2015)

fprintf([' \n *** Latency Difference: ' num2str(lat_in_ms) 'ms! *** \n'])

% Plot warping paths compared to diagonal
plotWP(i_QuerIntX, i_RefCorY, plot_GA_QuerIntX, plot_GA_RefCorY, analysis, timewindow, lat_in_ms, dtwpath);

%%             PERMUTATION TEST ON AREA MEASURE (have this for median(distance), i.e. lat_in_ms measure, too)
%       GET OBSERVED AREA
DIAG = 1:max(i_QuerIntX); % give me the diagonal
areaDIAG = trapz(DIAG); % area under diagonal
areaWP = trapz(i_QuerIntX, i_RefCorY);
area = (areaDIAG - areaWP) / areaDIAG; % area between diagonal and WP (whether this value is positive or negative informs about temporal order)

%       GET PERMUTED AREA
orig_QuerIntX = SSERP_Int;
orig_RefCorY = SSERP_Cor;

% Initialise
perm_x = cell(1,permutations); perm_y = cell(1, permutations); perm_area = cell(1, permutations);
% Update timepoints
timepoints = length(dataindices);
% Initialize permutation variables
for perms = 1:permutations
    perm_ERP1(1:timepoints, 1:subjects) = 0;
    perm_ERP2(1:timepoints, 1:subjects) = 0;
    for s = 1:subjects
        coin = randi(2, subjects, 1);  % flip the coin
        if coin(s) == 1
            perm_ERP1(:, s) = orig_RefCorY(:, s);  % assign based on coin
            perm_ERP2(:, s) = orig_QuerIntX(:, s);
        elseif coin(s) == 2
            perm_ERP1(:, s) = orig_QuerIntX(:, s);
            perm_ERP2(:, s) = orig_RefCorY(:, s);
        end
    end
    % Get Grand Averages
    perm_GA1 = mean(perm_ERP1, 2);
    perm_GA2 = mean(perm_ERP2, 2);
    % z-score permutation GAs
    perm_zGA1 = zscore(perm_GA1);
    perm_zGA2 = zscore(perm_GA2);
    
    % DTW
    [~, perm_x{perms}, perm_y{perms}] = dtw(perm_zGA1, perm_zGA2);
    
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
title({analysis, ['Latency diff.: ' num2str(lat_in_ms) 'ms - ', plot_p]}, 'FontSize', 20)
saveas(gcf, [dtwpath analysis ' -- PermDist.jpg'])
hold off
close all

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
title({analysis, ['Absolute Areas - Latency diff.: ' num2str(lat_in_ms) 'ms - ', plot_p]}, 'FontSize', 20)
saveas(gcf, [dtwpath analysis ' -- Absolute Areas -- PermDist.jpg'])
hold off
close all

%   Save permutation-variables to file
save([dtwpath 'DTW_PermTest_of_' analysis '.mat'], 'area', 'p', 'thresh1', 'thresh2', 'perm_x', 'perm_y', 'permdist_areas')
clear area p thresh
end

%% PLOT WARPING PATH (with non-standardized GAs on Axes)!
function plotWP(i_QuerIntX, i_RefCorY, plot_GA_QuerIntX, plot_GA_RefCorY, analysis, timewindow, lat_in_ms, dtwpath)
%% Plot warping path compared to diagonal
figure;
set(gcf,'Position', [ 0        0        1280         907]);
% WP
subplot(4,4, [2:4 6:8 10:12])
wp_fontsize = 17;
% Set some axes stuff
ax1 = gca;
ax1.XLim(2) = i_QuerIntX(end); ax1.YLim(2) = i_RefCorY(end);
tickstep = ax1.XTick(2)*2; % one step in ticks in milliseconds
ticks = timewindow(1)-2:tickstep:timewindow(2); % round down timewindow (1) from 302 to 300 
xmax = ax1.XLim(2);
hold on
% plot main diag
WPLineWidth = 3;
diag = plot(1:xmax, 1:xmax, 'r');
diag.Color = '#ff796c'; % salmon
diag.LineWidth = WPLineWidth;
% plot WP, 
p = plot(i_QuerIntX, i_RefCorY);
p.Color = '#8ab8fe'; % carolina blue
p.LineWidth = WPLineWidth;
% fill area between WP & main diagonal
f = fill(i_QuerIntX, i_RefCorY, rgb('pale green')); % pale green
f.FaceAlpha = 0.5;
f.EdgeAlpha = 0; % make this invisible
% this needs to be after plotting WP
ax1.XTickLabel = []; ax1.YTickLabel = []; % no tick labels for WP (but for GAs, below)

% Plot GA1 (reference / correct) on y-axis
subplot(4,4, [1 5 9])
pQ = plot(plot_GA_RefCorY);
pQ.Color = 'k'; pQ.LineWidth = 3;
ylabel('\muV')
camroll(90)
ax3 = gca;
ax3.XAxisLocation = 'top';

% Plot GA2 (query / intrusion) on x-axis
subplot(4,4, 14:16)
pR = plot(plot_GA_QuerIntX);
pR.Color = 'r'; pR.LineWidth = 3;
ylabel('\muV')
ax2 = gca;


ylimmax = max(max([plot_GA_QuerIntX, plot_GA_RefCorY])) + .5;
ax2.YLim(2) = ylimmax; ax3.YLim(2) = ax2.YLim(2);
ax2.XTick = ax1.XTick; ax3.XTick = ax2.XTick;
ax2.XLim(2) = length(plot_GA_QuerIntX); ax3.XLim(2) = ax2.XLim(2);
ax2.XTickLabel = ticks; ax3.XTickLabel = ticks; % set ticks dynamic depending on timewindow of ERPs
ax2.FontSize = wp_fontsize; ax3.FontSize = wp_fontsize;
% add some text in plots [needs to be different for N2pc because we are making it dependent from axes-limits]
if strcmp(analysis, 'N2pc')
    t1 = text(ax3, ax3.XLim(2) + 4, ylimmax - 0.5, 'Correct');
else
    t1 = text(ax3, ax3.XLim(2) + 10, ylimmax - 0.8, 'Correct');
end
if strcmp(analysis, 'N2pc')
    t2 = text(ax2, ax2.XLim(2) - 30, ylimmax + 0.28, 'Intrusion');
else
    t2 = text(ax2, ax2.XLim(2) - 60, ylimmax + 0.8, 'Intrusion');
end
t1.FontSize = 20; t2.FontSize = t1.FontSize;
t1.FontAngle = 'italic'; t2.FontAngle = t1.FontAngle;
t1.FontWeight = 'bold'; t2.FontWeight = t1.FontWeight;
hold off

% set superlabels
hold on
if strcmp(analysis, 'P3 at Fz')
    [~, h1] = suplabel(['Fz - Latency difference: ' num2str(lat_in_ms) 'ms'], 't', [.2 .11 .85 .85]);
else
    [~, h1] = suplabel([analysis ' - Latency difference: ' num2str(lat_in_ms) 'ms'], 't', [.2 .11 .85 .85]);
end
h1.FontSize = 22;
[~, h2] = suplabel('Time (ms)', 'x', [.08 .10 .85 .85]);
h2.FontSize = 20; 
[~, h3] = suplabel('Time (ms)', 'y');
h3.FontSize = h2.FontSize;
hold off

% Save it
saveas(gcf, [dtwpath analysis ' -- WarpPath.jpg'])
close all
end