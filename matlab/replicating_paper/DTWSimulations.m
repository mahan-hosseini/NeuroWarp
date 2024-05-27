%% FUNCTION DTWSimulations
%  Demonstrate validity of DTW approach for present dataset

function DTWSimulations(filepath, dtwpath)

% 00) convert strings to char
filepath = char(filepath);
dtwpath = char(dtwpath);

% 0A) Load Grand Average N2pc of Intrusions
loadstr = load([filepath 'LPF_ERP_structs.mat'], 'LPF_GA_N2pc');
TrueGA = loadstr.LPF_GA_N2pc.Intrusion;
CorrectGA = loadstr.LPF_GA_N2pc.Correct;
clear loadstr

% 0B) Global Variables
datapt_in_ms = 2;                                  % at 500Hz, 1 datapoint means 2ms
time = -100:datapt_in_ms:798;                      % our ERPs are sampled at 500Hz from -100 to 800ms
latencyshift = 50;                                 % shift of 50 ms
latinpts = latencyshift/datapt_in_ms;              % ShiftedGA has to be shifted by 25 datapoints to have 50ms shift in time
frames = length(TrueGA);                           % for noise function: number of datapts in signal
epochs = 1;                                        % for noise function: number of trials
srate = 500;                                       % for noise function: sampling rate
numruns = 25000;                                   % how often each noise level is run
ERPLineWidth = 1.5;

% 0C) Prepare plotpath to save plots
plotpath = [dtwpath 'DTW Simulations\'];
if ~isfolder(plotpath)
    mkdir(plotpath)
end

noisescaleranges = [0:0.05:0.95];                            

% only run simulations if we haven't before (if you want to run them again, delete the file)
if ~isfile([plotpath 'DTW Simulations Variables.mat']) % i know it's confusing that this is called plotpath but it's where we save everything 
    % noise scales: scaling of the noise added to the shifted N2pc
    noise_scales = noisescaleranges; 
    
    % store latency estimates
    latencyestimates = zeros(1,length(noise_scales));    % initialise latencyestimates
    latencyestimates_stds = zeros(1,length(noise_scales));% initialise std vals of lat estimates
    SNRs          = zeros(1,length(noise_scales));       % initialise SNRs
    
    for n_s = 1:length(noise_scales) % scalars used to regulate amplitude of noise added to shifted N2pc
        
        disp(['Noise Level: ' num2str(noise_scales(n_s))]) % some info while running
        
        these_SNRs = zeros(numruns,1); % repeat each noise level numruns times to decrease variability of randomness due to generation of noisy signal
        these_latestimates = zeros(numruns,1); % same for latency estimate
        for n = 1:numruns
            
            % 1) Add a fixed latency shift to GA N2pc Intrusion
            ShiftedGA = zeros(450,1);
            ShiftedGA(latinpts+1:end) = TrueGA(1:end-latinpts);
            
            % 2) Generate a Noise timeseries with Human Noise Code & add it to shifted GA N2pc Intrusion
            HumanNoise = noise(frames, epochs, srate);
            if size(HumanNoise,1) == 1
                HumanNoise = HumanNoise';
            end
            this_noise_scalar = noise_scales(n_s);
            ShiftedGA = HumanNoise * this_noise_scalar + ShiftedGA; % scale the noise we generate by the scalar & then just add it to the Shifted GA N2pc
            
            % 2b) A plot of true observed & noisy latency shifted GA N2pcs
            if n == 1
                if n_s == 1
                    ERPFig = figure;
                    set(gcf,'Position', [ 0        0        1280         907]);
                end
                if length(noise_scales) > 10
                    subplot(4,5,n_s)
                else
                    subplot(2,5,n_s)
                end
                TrueN2pcLine = plot(time, TrueGA, 'Color', 'r', 'LineWidth', ERPLineWidth);
                hold on
                ShiftedN2pcLine = plot(time, ShiftedGA, 'Color', '#2c6fbb', 'LineWidth', ERPLineWidth);
                ax1 = gca;
                ax1.XLim = [-100 800]; ax1.FontSize = 12;
                title(['Noisescale: ' num2str(this_noise_scalar)])
                hold off
            end
            
            % 3) DTW between shifted & true observed GA N2pc of Intrusions
            lat_in_ms = runDTW_local(ShiftedGA, TrueGA);
            
            % 4) Show how "off" DTW latency estimate is
            these_latestimates(n) = lat_in_ms;
            
            % 5) Calculate SNR using rootmeansquared(GA-timeseries(200:400)) / RMS(GA-timeseries(-100:100))
            % => RMS - How much am I away from zero in this timewindow
            noisewindow = [-50 100]; signalwindow = [200 400];                  % these are in ms - noise window starts at -50 and not -100 because -100:-50 is all zeros for shifted, which distorts comparison if included
            noiseindices = find(time>=noisewindow(1) & time<=noisewindow(2));    % find indices of our noise timewindow in datapoints
            signalindices = find(time>=signalwindow(1) & time<=signalwindow(2)); % same for signal
            these_SNRs(n) = rms(ShiftedGA(signalindices)) / rms(ShiftedGA(noiseindices)); % use noise & signal indices to compute root mean squared of those windows & divide them to get SNR
        end
        latencyestimates(n_s) = mean(these_latestimates);
        latencyestimates_stds(n_s) = std(these_latestimates);
        %         latencyerrors(n_s) = mean(these_latencyerrors);
        %         latencyerrors_SEs(n_s) = (std(these_latencyerrors) / sqrt(numruns));
        SNRs(n_s) = mean(these_SNRs);
    end
    all_latencyestimates = latencyestimates;
    all_latencyestimates_stds = latencyestimates_stds;
    all_SNRs = SNRs;

% 6) Make ERP Figure pretty & save it (we just do it whenever running,
%    because it's just showing example timeseries so for a given noise 
%    range we can re-use this figure)
figure(ERPFig)
[~, h1] = suplabel('Time (ms)', 'x', [.08 .1 .84 .84]);
[~, h2] = suplabel('Amplitude (\muV)', 'y', [.12 .08 .84 .84]);
[~, h3] = suplabel('Original (red) & noisy, time shifted (blue) Intrusion ERPs', 't', [.08 .13 .84 .84]);
h1.FontSize = 25; h2.FontSize = h1.FontSize; h3.FontSize = h1.FontSize;
saveas(ERPFig, [plotpath 'ERP Figure.png'])

% For plots - compute true observed SNR value of intrusion & correct N2pcs
TrueSNR = rms(TrueGA(signalindices)) / rms(TrueGA(noiseindices));
CorSNR = rms(CorrectGA(signalindices)) / rms(CorrectGA(noiseindices));

% save files for later overview plots (so we dont have to re-run
% simulations whenever we want to change something with plots)
save([plotpath 'DTW Simulations Variables'], 'all_latencyestimates', 'all_latencyestimates_stds', 'all_SNRs', 'TrueSNR', 'CorSNR'); % use plotpath even though we save vars (I'm a rebel what can you do)

else
    disp('there is a variable file - only plotting plots based on it') 
    % load variables 
    load([plotpath 'DTW Simulations Variables'], 'all_latencyestimates', 'all_latencyestimates_stds', 'all_SNRs', 'TrueSNR', 'CorSNR');
    
end

% plot main results
% Preparation - Get the correct variables for this analysis!
% A - The ones we saved
latencyestimates = all_latencyestimates;
latencyestimates_stds = all_latencyestimates_stds;
SNRs = all_SNRs;
% B - Noise Scales - scaling of the noise added to the shifted N2pc
noise_scales = noisescaleranges; 

% 8) Plot SNR as a function of Noise Scale
SNRScaleFig = figure;
SNRScaleColor = rgb('pine green');
plot(noise_scales, SNRs, 'Color', SNRScaleColor, 'LineWidth', 3);
if TrueSNR >= min(SNRs) && TrueSNR <= max(SNRs) % if SNR of original Intrusion N2pc is in current range, indicate its value
    SNRline = hline(TrueSNR);
    CorSNRline = hline(CorSNR);
    SNRline.Color = SNRScaleColor; CorSNRline.Color = SNRScaleColor;
    SNRline.LineWidth = 1; CorSNRline.LineWidth = 1;
    SNRline.LineStyle = '--'; CorSNRline.LineStyle = '--';
end
ax = gca;
ax.FontSize = 18;
ax.Title.String = 'SNR as a function of noise'; ax.Title.FontSize = 20;
ax.XLabel.String = 'Noise Level (scalar used on noise)';
ax.YLabel.String = 'Signal to Noise Ratio (SNR)';
ax.XLim = [min(noise_scales) max(noise_scales)];
ax.YLim = [0 20]; % to make sure that there is room to indicate correct N2pc's SNR
saveas(SNRScaleFig, [plotpath 'SNR by Noise Figure.png'])

% preparation 2 - create vars that have +- std around latency estimates
latencyestimates_SE = latencyestimates_stds / sqrt(numruns);
latency_plusSE = latencyestimates + latencyestimates_SE;
latency_minusSE = latencyestimates - latencyestimates_SE;

% 9) Plot DTW accuracy as a function of Noise Scale
AccScaleFig = figure;
fillX = [noise_scales, fliplr(noise_scales)];
fillY = [latency_minusSE, fliplr(latency_plusSE)];
AccScaleColor = rgb('royal purple');
patch(fillX, fillY, AccScaleColor, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
hold on
plot(noise_scales, latencyestimates, 'Color', AccScaleColor, 'LineWidth', 3)
hold off
ax = gca;
ax.FontSize = 18;
ax.Title.String = 'Latencies as a function of noise'; ax.Title.FontSize = 20;
ax.XLabel.String = 'Noise Level (scalar used on noise)';
ax.YLabel.String = 'Latency Estimate (ms)';
saveas(AccScaleFig, [plotpath 'DTW Accuracy by Noise Figure.png'])

% 10) Plot DTW accuracy as a function of SNR
AccSNRFig = figure;
fillX = [SNRs, fliplr(SNRs)];
fillY = [latency_minusSE, fliplr(latency_plusSE)];
AccSNRColor = rgb('auburn');
patch(fillX, fillY, AccSNRColor, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
hold on
plot(SNRs, latencyestimates, 'Color', AccSNRColor, 'LineWidth', 3)
set(gca, 'XDir', 'reverse')
ax = gca;
ax.XLim = [min(SNRs) max(SNRs)];
ax.FontSize = 18;
ax.Title.String = 'Latencies as a function of SNR'; ax.Title.FontSize = 20;
ax.XLabel.String = 'Signal to Noise Ratio (SNR)';
ax.YLabel.String = 'Latency Estimate (ms)';
if TrueSNR >= min(SNRs) && TrueSNR <= max(SNRs) % if SNR of original Intrusion N2pc is in current range, indicate its value
    SNRline = vline(TrueSNR);
    SNRline.Color = AccSNRColor; 
    SNRline.LineWidth = 1; 
    SNRline.LineStyle = '--'; 
end
hold off
saveas(AccSNRFig, [plotpath 'DTW Accuracy by SNR.png'])
end

%% Run DTW locally between True observed GA N2pc & Shifted GA N2pc
function lat_in_ms = runDTW_local(ShiftedGA_queryX, TrueGA_refY)
%       ---- MAIN DTW COMPUTATION HAPPENS HERE ----
% Outputs IX & IY == indices of x & y coordinates of warping path!
%   ORIGINAL FORM: [dist, ix, iy] = dtw(x,y)
%   ==> QUERY/SHIFTED IS ON X-AXIS OF WP, SO FIRST INPUT
%   ==> REF/TRUE IS ON Y-AXIS OF WP, SO SECOND INPUT!
%   ==> WP of this will be below main diagonal if ShiftedGA (QUERY) is
%       LATER than TrueObserved N2pc (REFERENCE)
%   ==> I.e. AREA will be positive in this case because area = areaDIAG -
%       areaWP (and hence areaDIAG will be larger than areaWP if WP is below
%       DIAG [and ShiftedGA (Q) is indeed later than TrueGA(R)!]

% Only test the timewindow we have in N2pcP3 correlation code for the N2pc too
time = -100:2:798;
N2pctimewin = [150 398];
N2pctime = N2pctimewin(1):2:N2pctimewin(2); % because EEG @ 500 Hz
N2pcinds = find(time>=N2pctime(1) & time <=N2pctime(end));

ShiftedGA_queryX = ShiftedGA_queryX(N2pcinds); % extract time window of interest
TrueGA_refY = TrueGA_refY(N2pcinds);

ShiftedGA_queryX = zscore(ShiftedGA_queryX); % standardized before running DTW
TrueGA_refY = zscore(TrueGA_refY);

[~, i_ShiftedGA_queryX, i_TrueGA_refY] = dtw(ShiftedGA_queryX, TrueGA_refY);

% Latency in MS is distribution of point-wise distances (distance = ix - iy) / (sampling rate (Hz) / 1000)
distance = i_ShiftedGA_queryX - i_TrueGA_refY;
latency = median(distance);
lat_in_ms = latency / .5; % 500Hz
% this is as proposed in Zoumpoulaki et al. (2015)
end

%% Noise from the Human EEG Frequency Spectrum
% From: Yeung N, Bogacz R, Holroyd CB, Cohen JD. Detection of synchronized oscillations in the electroencephalogram: an evaluation of methods. Psychophysiology. 2004;41(6):822–832
function signal = noise(frames, epochs, srate)

% Function generates noise with the power spectrum of human EEG
% Inputs:
%  frames - number of signal frames per each trial
%  epochs - number of simulated trials
%  srate - sampling rate of simulated signal
% Output:
%  signal - simulated EEG signal; vector: 1 by frames*epochs containing concatenated trials
% Implemented by: Rafal Bogacz and Nick Yeung, Princeton Univesity, December 2002

load meanpower
sumsig = 20;	%number of sinusoids from which each simulated signal is composed of

signal = zeros (1, epochs * frames);
for trial = 1:epochs
    freq=0;
    range = [(trial-1)*frames+1:trial*frames];
    for i = 1:sumsig
        freq = freq + (4*rand(1));
        freqamp = meanpower(min (ceil(freq), 125)) / meanpower(1);
        phase = rand(1)*2*pi;
        signal (range) = signal (range) + sin ([1:frames]/srate*2*pi*freq + phase) * freqamp;
    end
end
end
