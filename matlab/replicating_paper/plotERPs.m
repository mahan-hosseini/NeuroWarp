function plotERPs(P3s, N2pcs, plotpath, LPF, plotSE)
%%  Initialize
% timepoints var
tps = 1:size(P3s(1).Correct,1);
%timelims = [-100 798];
subjects = 1:length(P3s);
P3chans = {'Fz', 'Cz', 'Pz'};

% convert strings to char if needed
plotpath = char(plotpath);

% experiment 2?
if length(subjects) < 15
    exp2 = 1;
else
    exp2 = 0;
end

% fix plotpath based on LPF & exp2
if LPF
    plotpath = [plotpath 'LPF_ERPs/'];
    if exp2
        plotpath = [plotpath 'Exp2/'];
    end
else
    plotpath = [plotpath 'Original_ERPs/'];
    if exp2
        plotpath = [plotpath 'Exp2/'];
    end
end
% check if plotpath exists, if not create it
if ~isfolder(plotpath)
    mkdir(plotpath)
end

% axis-settings
% SS_FSize = 17;
SS_FSize = 13; % laptop

% P3s
for c = 1:length(P3chans)
    %   SS P3s
    figure;
    ymin = 0; ymax = 0;
    for s = 1:length(subjects)
        this_ymin(s) = min(vertcat(P3s(s).Correct(:,c), P3s(s).Intrusion(:,c)));
        this_ymax(s) = max(vertcat(P3s(s).Correct(:,c), P3s(s).Intrusion(:,c)));
    end
    ymin = min(this_ymin); ymax = max(this_ymax);
    set(gcf,'Position', [ 0        0        1280         907]);
    for s = 1:length(subjects)
        subplot(5, 5, s)
        p = plot(tps, P3s(s).Correct(:,c));
        p.LineWidth = 1.3; p.Color = 'k';
        hold on
        title(['SUBJ #' num2str(subjects(s))])
        p = plot(tps, P3s(s).Intrusion(:,c));
        p.LineWidth = 1.3; p.Color = 'r';
        xlim([1 tps(end)])
        xticks([1 50 100 150 200 250 300 350 400 450])
        xticklabels(-100:100:800)
        xtickangle(45)
        ylim([ymin ymax])
        h = hline(0, 'k:');
        h.LineWidth = 1.5;
        set(gca,'FontSize',SS_FSize)
        hold off
    end
    hold on
    [~,h1] = suplabel(P3chans{c},'t');
    %     set(h1, 'FontSize', 35);
    set(h1, 'FontSize', 20);
    hold off
    if LPF
        saveas(gcf, [plotpath 'LPF_SS_ ' P3chans{c} '.jpg'])
    else
        saveas(gcf, [plotpath 'SS_ ' P3chans{c} '.jpg'])
    end
    
    %   GA P3s
    figure;
    set(gcf,'Position', [ 0        0        1280         907]);
    hold on
    tmp1 = cat(3,P3s.Correct);
    SSERPs = squeeze(tmp1(:, c, :));
    plotSEtoo(SSERPs,'k', plotSE)
    hold off
    hold on
    tmp1 = cat(3,P3s.Intrusion);
    SSERPs = squeeze(tmp1(:, c, :));
    plotSEtoo(SSERPs,'r', plotSE)
    xlim([1 tps(end)])
    xticks([1 50 100 150 200 250 300 350 400 450])
    xticklabels(-100:100:800)
    ylim([-5 10])
    h = hline(0, 'k:');
    h.LineWidth = 2.5;
        set(gca,'FontSize', 40)
%     set(gca, 'FontSize', 25);
    [~,h1] = suplabel(P3chans{c},'t');
        set(h1, 'FontSize', 45)
%     set(h1, 'FontSize', 30)
    hold off
    if LPF
        saveas(gcf, [plotpath 'LPF_GA_ ' P3chans{c} '.jpg'])
    else
        saveas(gcf, [plotpath 'GA_ ' P3chans{c} '.jpg'])
    end
end

%% N2pcs [NOTE - WE HAVE DIFFERENT XLIMS & TICKS FOR N2pc (only -100 to 500)!]
%   SS N2pcs
figure;
set(gcf,'Position', [ 0        0        1280         907]);
ymin = 0; ymax = 0; this_ymin = 0; this_ymax = 0;
for s = 1:length(subjects)
    this_ymin(s) = min([min(N2pcs(s).Correct) min(N2pcs(s).Intrusion)]);
    this_ymax(s) = max([max(N2pcs(s).Correct) max(N2pcs(s).Intrusion)]);
end
ymin = min(this_ymin); ymax = max(this_ymax);
for s = 1:length(subjects)
    subplot(5, 5, s)
    p = plot(tps, N2pcs(s).Correct);
    p.LineWidth = 1.3; p.Color = 'k';
    hold on
    title(['SUBJ #' num2str(subjects(s))])
    p = plot(tps, N2pcs(s).Intrusion);
    p.LineWidth = 1.3; p.Color = 'r';
    xlim([1 300])
    xticks([1 50 100 150 200 250 300])
    xticklabels(-100:100:500)
    xtickangle(45)
    ylim([ymin ymax])
    h = hline(0, 'k:');
    h.LineWidth = 1.5;
    set(gca,'FontSize',SS_FSize)
    hold off
end
hold on
[~,h1] = suplabel('N2pc','t');
%     set(h1, 'FontSize', 35);
set(h1, 'FontSize', 20);
hold off
if LPF
    saveas(gcf, [plotpath 'LPF_SS_N2pc.jpg'])
else
    saveas(gcf, [plotpath 'SS_N2pc.jpg'])
end

%   GA N2pcs
% If Exp 1, do one GA N2pc plot for part A, one for part B and one combined
if ~exp2 
    % Experiment 1A
    Exp1A = N2pcs(1:12);
    Corrects1A = cat(3,Exp1A.Correct);
    Corrects1A = squeeze(Corrects1A);
    Intrusions1A = cat(3,Exp1A.Intrusion);
    Intrusions1A = squeeze(Intrusions1A);
    
    % plot
    figure;
    set(gcf,'Position', [ 0        0        1280         907]);
    hold on
    plotSEtoo(Corrects1A,'k', plotSE)
    hold off
    hold on
    plotSEtoo(Intrusions1A,'r', plotSE)
    hold on
    xlim([1 300])
    xticks([1 50 100 150 200 250 300])
    xticklabels(-100:100:500)
    ylim([-5 1])
    yticks(-5:1:1)
    h = hline(0, 'k:');
    h.LineWidth = 2.5;
    set(gca, 'FontSize', 40)
    [~,h1] = suplabel('N2pc - Exp 1A','t');
    set(h1, 'FontSize', 45)
    hold off
    if LPF
        saveas(gcf, [plotpath 'LPF_GA_N2pc_1A.jpg'])
    else
        saveas(gcf, [plotpath 'GA_N2pc_1A.jpg'])
    end
    
    % Experiment 1B
    Exp1B = N2pcs(13:end);
    Corrects1B = cat(3,Exp1B.Correct);
    Corrects1B = squeeze(Corrects1B);
    Intrusions1B = cat(3,Exp1B.Intrusion);
    Intrusions1B = squeeze(Intrusions1B);
    
    % plot
    figure;
    set(gcf,'Position', [ 0        0        1280         907]);
    hold on
    plotSEtoo(Corrects1B,'k', plotSE)
    hold off
    hold on
    plotSEtoo(Intrusions1B,'r', plotSE)
    hold off
    xlim([1 300])
    xticks([1 50 100 150 200 250 300])
    xticklabels(-100:100:500)
    ylim([-5 1])
    yticks(-5:1:1)
    h = hline(0, 'k:');
    h.LineWidth = 2.5;
    set(gca, 'FontSize', 40)
    [~,h1] = suplabel('N2pc - Exp 1B','t');
    set(h1, 'FontSize', 45)
    hold off
    if LPF
        saveas(gcf, [plotpath 'LPF_GA_N2pc_1B.jpg'])
    else
        saveas(gcf, [plotpath 'GA_N2pc_1B.jpg'])
    end
end

figure;  % this is the (experiments-combined) figure we show in the paper
set(gcf,'Position', [ 0        0        1280         907]);
tmp1 = cat(3,N2pcs.Correct);
SSERPs = squeeze(tmp1);
SSERPs = SSERPs';
hold on
plotSEtoo(SSERPs,'k', plotSE)
hold off
hold on
tmp1 = cat(3,N2pcs.Intrusion);
SSERPs = squeeze(tmp1);
plotSEtoo(SSERPs,'r', plotSE)
xlim([1 300])
xticks([1 50 100 150 200 250 300])
xticklabels(-100:100:500)
ylim([-5 1])
yticks(-5:1:1)
h = hline(0, 'k:');
h.LineWidth = 2.5;
set(gca,'FontSize', 40)
[~,h1] = suplabel('N2pc','t');
    set(h1, 'FontSize', 45)
hold off
if LPF
    saveas(gcf, [plotpath 'LPF_GA_N2pc.jpg'])
else
    saveas(gcf, [plotpath 'GA_N2pc.jpg'])
end
end

%% plotSEtoo (plot SE around mean for Grand Averages as well)
function plotSEtoo(SS_ERP, color, plotSE) %takes single-subject ERPs (subj x tps) and color

% if matrix is the wrong way around (tps x trialnums instead of trials x
% tps) transpose it!
if size(SS_ERP,1) > size(SS_ERP,2) 
    SS_ERP = SS_ERP';
end

trialnums = size(SS_ERP,1);
timepoints = size(SS_ERP,2);
sd = std(SS_ERP);
SE = sd / sqrt(trialnums);
s1 = mean(SS_ERP,1)-SE;
s2 = mean(SS_ERP,1)+SE;
GA = mean(SS_ERP,1);

% hold on
p1 = plot(1:timepoints,GA);
p1.LineWidth = 4; p1.Color = color;
fillX = [1:timepoints, ...
    fliplr(1:timepoints)];
fillY = [s1, fliplr(s2)];
if plotSE
    f1 = fill(fillX, fillY, color, 'HandleVisibility', 'off'); %MAHAN: handle_vis = off so it won't be included in legend
    f1.FaceAlpha = 0.15; f1.EdgeColor = color;
    f1.EdgeAlpha = 0.15; f1.LineStyle = ':';
end
hold off
end
