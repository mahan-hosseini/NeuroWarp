% How to define latency correlations between N2pc and P3
%   1) Get Single Subject (or group) ERPs
%   2) Get N bootstrap ERPs (w. same bootstrapping subject-distributions 
%       for N2pc as well as P3)
%   3) DTW of 2) with 1) using AREA measure!
%   4) Gives you two numbers for each bootstrap sample -> DTW
%       result of N2pc & P3 latency (lag)
%   5) Correlate these with Pearson and Spearman (latter preferred for
%       presentations)

function correlateP3andN2pc(filepath, P3timewin, N2pctimewin, num_boots, outlier)
%%                          PREPARATION

% strings to char
filepath = char(filepath);

% load 
load([filepath 'LPF_ERP_structs.mat'])

% since we use smoothed ERPs, rename them for DTW functions below
P3s = LPF_P3s;
N2pcs = LPF_N2pcs;

conditions = {'Correct', 'Intrusion'};
P3str = {'Fz', 'Cz', 'Pz'}; %again, 1 = Fz, 2 = Cz, 3 = Pz
% prepare path for plots
plotpath = [filepath 'Plots/Latency_correlations/'];
if ~exist(plotpath,'dir')
    mkdir(plotpath)
end

full_boots = num_boots;
savediaryto = [plotpath 'CorrelateAnalysis_Info.txt'];
diary(savediaryto)

% for saving marginals & stats to file
N2pc_marginals = struct('Correct', [], 'Intrusion', []);
P3_marginals = struct('Correct', [], 'Intrusion', []);
Stats = struct('Correct', [], 'Intrusion', []);

for c = 1:length(conditions)
    
    %%          THIS BIT MAKES IT SO THAT WE ONLY ANALYISE PZ
    % Not the prettiest approach because this function originally looped
    % over channels 
    P3chan = 3; % ONLY PZ IN THIS FUNCTION !!!
    p = 3;      % ONLY PZ IN THIS FUNCTION !!!
    
    %%                  RUN THINGS
    % loop variables are inputs in local function
    condition = conditions{c};
    % reset number of bootstrap (changed in code below) to be number of bootstrap you wanted in outer function
    num_boots = full_boots;
    % time-window you want for each component
    if ~isequal(size(P3timewin), [1 2]) % allow P3 time window to be different for different channels (has to be chan x 2 size!)
        P3time = P3timewin(p,1):2:P3timewin(p,2);
    else
        P3time = P3timewin(1):2:P3timewin(2);
    end
    N2pctime = N2pctimewin(1):2:N2pctimewin(2);
    % run DTW & bootstrap
    [P3_dtw_area, N2pc_dtw_area, all_boots_P3s, all_boots_N2pcs] = runbootsDTW(P3s,N2pcs,condition,P3chan,P3time,N2pctime,num_boots);
    
    %% OUTLIER REJECTION? 
    if outlier
        % exclude outliers if +-X SD

        P3_std = std(P3_dtw_area);
        N2pc_std = std(N2pc_dtw_area);
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
        P3pos_outlier_idx = find(P3_dtw_area > (mean(P3_dtw_area) + (std_thresh*P3_std)));
        if ~isempty(P3pos_outlier_idx)
            outlier_idx = vertcat(outlier_idx, P3pos_outlier_idx);
        end
        P3neg_outlier_idx = find(P3_dtw_area < (mean(P3_dtw_area) - (std_thresh*P3_std)));
        if ~isempty(P3neg_outlier_idx)
            outlier_idx = vertcat(outlier_idx, P3neg_outlier_idx);
        end
        N2pcpos_outlier_idx = find(N2pc_dtw_area > (mean(N2pc_dtw_area) + (std_thresh*N2pc_std)));
        if ~isempty(N2pcpos_outlier_idx)
            outlier_idx = vertcat(outlier_idx, N2pcpos_outlier_idx);
        end
        N2pcneg_outlier_idx = find(N2pc_dtw_area < (mean(N2pc_dtw_area) - (std_thresh*N2pc_std)));
        if ~isempty(N2pcneg_outlier_idx)
            outlier_idx = vertcat(outlier_idx, N2pcneg_outlier_idx);
        end
        outlier_idx = unique(outlier_idx); % get rid of duplicates
        
        P3_dtw_area(outlier_idx) = []; % remove outlier indices
        N2pc_dtw_area(outlier_idx) = [];
        
    end
    
    %%         Save Var, Skew & Kurtosis Values of DTW Area Distributions
    Stats.(conditions{c}) = struct('P3Var', var(P3_dtw_area), 'P3Skew', skewness(P3_dtw_area), 'P3Kurt', kurtosis(P3_dtw_area), ...
        'N2pcVar', var(N2pc_dtw_area), 'N2pcSkew', skewness(N2pc_dtw_area), 'N2pcKurt', kurtosis(N2pc_dtw_area));
    if outlier
        save([filepath num2str(std_thresh) 'SD OutlierRejected Marginal Stats.mat'], 'Stats') % save marginal-stats to file
    else
        save([filepath 'Marginal Stats.mat'], 'Stats')
    end
    
    %%              Standardise DTW Area Distributions 
    % ==> So line of linear fit in scatter & marginal dists reflect
    %   correlation values better (since correlations have standardisation
    %   internally)
    x = P3_dtw_area;
    y = N2pc_dtw_area;
    
    % standardize x & y
    x = zscore(x);
    y = zscore(y);   
    
    z_x= x; % for saving standardised and original marginals later
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
    xlimits = [-8 8]; ylimits = [-5 5]; % when zscoring marginals
    ylimmax = 800;
    %% first, scatterplot
    subplot(4,4,[2:4 6:8 10:12])
    scatter(x,y, [], scatter_color);
    hold on
    ax = gca;
    ax.FontSize = ax_fontsize;
    ax.YLim = ylimits; % when zscoring marginals
    ax.XLim = xlimits;
    P3_xlimits = ax.XLim;   % use the scatterplot's xlimits for P3's marginal xlimits
    N2pc_xlimits = ax.YLim; % use the scatterplot's Y (!) limits for N2pc's marginal X (!) limits
    ax.XTickLabels = [];
    ax.YTickLabels = [];
    title({'DTW area-diff correlations of P3 & N2pc.', ...
        ['For ' condition ' trials & P3 at ' P3str{P3chan} '.'], ...
        ['Pearson: r(' str_df ') = ' pears_str_r ': ' pears_str_p ' & Spearman: r(' str_df ') = ' ...
        spear_str_r ': ' spear_str_p]}, 'fontsize', 20)
    l = lsline;
    l.LineWidth = 3;
    l.Color = line_color;
    %% second, vN2pc marginal distribution
    % initialise axes
    subplot(4,4,[1 5 9])
    ax_N2pc = gca;
    ax_N2pc.FontSize = ax_fontsize;
    ax_N2pc.XLim = N2pc_xlimits;
    ax_N2pc.YLim(2) = ylimmax;
    camroll(90)
    hold on
    % reference lines on marginal & fix yticks
    reflines = hline([200 400 600], {'k', 'k', 'k'}); % reference lines for comparison
    for r = 1:length(reflines)
        reflines(r).LineWidth = 0.2;
        reflines(r).Color = rgb('teal');
        reflines(r).LineStyle = '--';
    end
    ax_N2pc.YTick = 0:200:800;
    ax_N2pc.YTickLabels = {'', '200', '400', '600', ''};
    ax_N2pc.YTickLabelRotation = 45;
    % histogram
    h_N2pc = histogram(y, 100);
    h_N2pc.FaceColor = scatter_color;
    h_N2pc.EdgeColor = line_color;
    set(ax_N2pc,'xaxislocation', 'top') % put xticklabels on top
%     % variance, skewness & kurtosis                                       % INSTEAD OF PUTTING THESE IN THE TITLE, I WILL USE THE VALUES SAVED TO FILE (above) AND CREATE A TABLE IF WANTED!
%     var_N2pc = var(y);
%     var_N2pc = round(var_N2pc,3);
%     skew_N2pc = skewness(y); % compute skewness
%     skew_N2pc = round(skew_N2pc,3);
%     kurt_N2pc = kurtosis(y); % compute kurtosis
%     kurt_N2pc = kurt_N2pc - 3; % standardize it so if kurt == normal, it's 0
%     kurt_N2pc = round(kurt_N2pc,3);
%     ax_N2pc.XLabel.String = {'N2pc Marginal Distribution', ['Var: ' num2str(var_N2pc) ' | Skew: ' num2str(skew_N2pc) ' | Kurt: ' num2str(kurt_N2pc)]};
    ax_N2pc.XLabel.String = 'N2pc Marginal Distribution';
    ax_N2pc.XLabel.FontSize = 17; % "Title" as XLabel because of axis-rotation
    ax_N2pc.XLabel.FontWeight = 'bold';
    hold off
    %% third, vP3 marginal distribution
    subplot(4,4,14:16)
    % initialise axes
    ax_P3 = gca;
    ax_P3.FontSize = ax_fontsize;
    ax_P3.XLim = P3_xlimits;
    ax_P3.YLim(2) = ylimmax;
    hold on
    % reference lines & yticks
    reflines = hline([200 400 600], {'k', 'k', 'k'}); % reference lines for comparison
    for r = 1:length(reflines)
        reflines(r).LineWidth = 0.2;
        reflines(r).Color = rgb('teal');
        reflines(r).LineStyle = '--';
    end
    ax_P3.YTick = 0:200:800;
    ax_P3.YTickLabels = {'', '200', '400', '600', ''};
    ax_P3.YAxisLocation = 'right';
    % histogram
    h_P3 = histogram(x, 100);
    h_P3.FaceColor = scatter_color;
    h_P3.EdgeColor = line_color;
    % variance, skewness & kurtosis
%     var_P3 = var(x);       % compute variance
%     var_P3 = round(var_P3,3);
%     skew_P3 = skewness(x); % compute skewness
%     skew_P3 = round(skew_P3,3);
%     kurt_P3 = kurtosis(x); % compute kurtosis
%     kurt_P3 = kurt_P3 - 3; % standardize it so if kurt == normal, it's 0
%     kurt_P3 = round(kurt_P3,3);
%     title({'P3 Marginal Distribution', ['Var: ' num2str(var_P3) ' | Skew: ' num2str(skew_P3) ' | Kurt: ' num2str(kurt_P3)]}, 'FontSize', 15)
    title('P3 Marginal Distribution', 'FontSize', 15)
    hold off
    %% Reposition axes to remove space between subplots, add labels & save
    scatterPos = ax.Position;
    N2pcPos = ax_N2pc.Position;
    P3Pos = ax_P3.Position;
    N2pcPos(1) = scatterPos(1) - N2pcPos(3);
    P3Pos(2) = scatterPos(2) - P3Pos(4);
    ax_N2pc.Position = N2pcPos;
    ax_P3.Position = P3Pos;
%    EXCLUDED -  % super labels
%     hX = suplabel('(DTW-areas between) True observed & surrogate P3s', 'x', [.21 .18 .84 .84]);
%     hX.FontSize = 22;
%     hY = suplabel('(DTW-areas between) True observed & surrogate N2pc', 'y', [.14 .22 .84 .84]);
%     hY.FontSize = 22;
    % save it
    if ~exist(plotpath, 'dir')
        mkdir(plotpath)
    end
    if outlier
        saveas(gcf,[plotpath num2str(std_thresh) 'SD_OutlierRemoved_' conditions{c} '_' P3str{p} '_&_N2pc_NoOut_LatCorr.jpg'])
    else
        saveas(gcf,[plotpath conditions{c} '_' P3str{p} '_&_N2pc_NoOut_LatCorr.jpg'])
    end
    hold off
    % save marginals to struct
    N2pc_marginals.(conditions{c}) = y;
    P3_marginals.(conditions{c}) = x;
    zN2pc_marginals.(conditions{c}) = z_y;
    zP3_marginals.(conditions{c}) = z_x;
end
if outlier
    save([filepath num2str(std_thresh) 'SD OutlierRejected Marginal Distributions.mat'], 'N2pc_marginals', 'P3_marginals', 'zN2pc_marginals', 'zP3_marginals') % save marginals to file
else
    save([filepath 'Marginal Distributions.mat'], 'N2pc_marginals', 'P3_marginals', 'zN2pc_marginals', 'zP3_marginals')
end
diary off
end

%%                      RUN THE DTW BOOTSTRAP STUFF
function [P3_dtw_area, N2pc_dtw_area, all_boots_P3s, all_boots_N2pcs] = ...
    runbootsDTW(P3s,N2pcs,condition,P3chan, P3time, N2pctime, num_boots)

% extract SS-components & GA from structures
[P3, N2pc, zP3, zN2pc, GA_P3, GA_N2pc, zGA_P3, zGA_N2pc] = extractSSGA(P3s, N2pcs, condition, P3chan, P3time, N2pctime);

% bootstrap manually using randi & run dtw between bootstrap GA & observed GA
num_subs = length(P3s);
num_bootstraps = num_boots;
P3_dtw_area(num_bootstraps,1) = 0; N2pc_dtw_area(num_bootstraps,1) = 0; % initialize stuff
% and loop
for n = 1:num_bootstraps
    
    disp(['boots num #' num2str(n)]) % info
    
    % get num_subj bootstrap samples (using randi == with duplicates)
    this_trap = randi(num_subs,num_subs,1); % max-val = num_subs, size = num_subs x 1
    boot_P3 = P3(:, this_trap);
    boot_N2pc = N2pc(:, this_trap);
    
    GA_boot_P3 = mean(boot_P3,2);
    GA_boot_N2pc = mean(boot_N2pc,2);
    
    zGA_boot_P3 = zscore(GA_boot_P3);
    zGA_boot_N2pc = zscore(GA_boot_N2pc);
    
    % DTW between GA_boot & GA_obs
    %   initialize
    ix = []; iy = []; distance = [];
    
    %       P3
    %   IX & IY ARE indices of x & y coordinates of warping path!
    [~, ix, iy] = dtw(zGA_boot_P3, zGA_P3);
    DIAG = 1:max(ix); % get the diagonal
    areaDIAG = trapz(DIAG); % area under diag
    areaWP = trapz(ix, iy); % area under WP
    P3_dtw_area(n) = (areaDIAG - areaWP) / areaDIAG; % to standardize, as in Zoumpoulaki et al 2015
    
    %   reset
    clear ix iy areaWP areaDIAG
    
    %       N2pc
    [~, ix, iy] = dtw(zGA_boot_N2pc, zGA_N2pc);
    DIAG = 1:max(ix); % get the diagonal
    areaDIAG = trapz(DIAG); % area under diag
    areaWP = trapz(ix, iy); % area under WP
    N2pc_dtw_area(n) = (areaDIAG - areaWP) / areaDIAG;
    
    %   Save all bootstraped grand averages for investigations
    all_boots_P3s{n} = zGA_boot_P3;
    all_boots_N2pcs{n} = zGA_boot_N2pc;
end
end

%%                EXTRACT SingleSubj GAs & zscored ones
function [P3, N2pc, zP3, zN2pc, GA_P3, GA_N2pc, zGA_P3, zGA_N2pc] = extractSSGA(P3s, N2pcs, condition, P3chan, P3time, N2pctime)
% narrow time-window of interest: P3 = 250 to 800 / N2pcs 150 to 400
time = -100:2:798;

P3inds = find(time>=P3time(1) & time<=P3time(end));
N2pcinds = find(time>=N2pctime(1) & time <=N2pctime(end));

tmp = cat(3,P3s.(condition));
tmp2 = tmp(:,P3chan,:);
full_P3 = squeeze(tmp2);

clear tmp tmp2
tmp = cat(3,N2pcs.(condition));
full_N2pc = squeeze(tmp);

P3 = []; N2pc = [];
num_subs = length(P3s);
for i = 1:num_subs
    P3(:, i) = full_P3(P3inds,i);
    N2pc(:, i) = full_N2pc(N2pcinds,i);
end

% need observed GAs, too
GA_P3 = mean(P3,2);
GA_N2pc = mean(N2pc,2);

% z-score P3 & N2pc
for i = 1:num_subs
    zP3(:, i) = zscore(P3(:, i));
    zN2pc(:, i) = zscore(N2pc(:, i));
end

% % and z-scored GAs
zGA_P3 = zscore(GA_P3);
zGA_N2pc = zscore(GA_N2pc);
end