%% Clear everything and establish where data is
clearvars;
clc;
if exist('E:\Documents\Recorded Data\EEG Feb 2015', 'dir') % location on desktop
    file_directory = 'E:\Documents\Recorded Data\EEG Feb 2015';
elseif exist('D:\Recorded Data', 'dir') % location on laptop
    file_directory = 'D:\Recorded Data';
elseif exist('C:\EEG Data\mit-data', 'dir')
    file_directory = 'C:\EEG Data\mit-data';
elseif exist('/Users/jacksonclee/Dropbox (MIT)/Projects/EEG_Frequencies/runScripts/Results', 'dir')
    file_directory = '/Users/jacksonclee/Dropbox (MIT)/Projects/EEG_Frequencies/runScripts/Results';
elseif exist('/Users/robertsonce/Dropbox (MIT)/Projects/EEG_Frequencies/runScripts/Results', 'dir')
    file_directory = '/Users/robertsonce/Dropbox (MIT)/Projects/EEG_Frequencies/runScripts/Results';
else
    error('please provide directory where file is stored');
end

% setFigProps

%% File Names
n{1} = 1;
analRunNum = 4;

%% Analysis Variables
discard_start = 0; % how much time should be cut off at beginning
occipitals = [20 26 27 28 29 30 31];
electrodes = occipitals;%[1:32]; %[27, 29, 64]; [27, 29, 30,28];
%electrodes = [1:32]; %[27, 29, 64]; [27, 29, 30,28];

%stimulation_frequencies = [28.33]; %, 21.25];
% intermodulation_frequencies = [21.6, 43.2];
% intermodulation_frequencies = 21.6;

% Dates & Day #: day2 - 3_11; day3 - 3_14; day4 - 3_15; day5 - 3_17; day6 -
% 3_21

fileNames = dir('../../runScripts/Results/stratus01*3_28*');
fileNames = {fileNames.name};
filenames{1}(1).name = char(fileNames(analRunNum));
legend = importdata('pilotlegends/stratus01.txt');

stimulation_frequencies = legend.data(analRunNum,1:2);
trialNum = legend.data(analRunNum,3);
trial_dur = legend.data(analRunNum,4); %JL  16; %JF

% add intermodulation freq as another freq of interest
% stimulation_frequencies(3) = (2*min(stimulation_frequencies))-max(stimulation_frequencies);

% group_freqs = cell(2, max([n{1}, n{2}]));
focus_electrode = 31; %27

%% Practice Analysis with one subject
fileID = fullfile(file_directory, filenames{1}(1).name);

%% load config params
[cfg_trldef, cfg_preproc] = configParams(fileID,discard_start,trial_dur,trialNum,electrodes)

%% preprocess data
all_data = ft_preprocessing(cfg_preproc);


%% Percept Definition
% Find all Percept Starts
% cfg_percdef.dataset = fileID;
% cfg_percdef.trialdef.eventtype = 'STATUS';
% cfg_percdef.trialfun = 'ft_trialfun_general';
% cfg_percdef.trialdef.prestim = 0;
% cfg_percdef.trialdef.poststim = 0;
% cfg_percdef.trialdef.eventvalue = 1:10;
% cfg_percdef = ft_definetrial(cfg_percdef);
% 
% 
% for iTrial = 1:trialNum
%     
%     trial_start = cfg_trldef.trl(iTrial, 1);
%     trial_end = cfg_trldef.trl(iTrial, 2);
%     
%     currPercepts = (cfg_percdef.trl(:, 1) >= trial_start) & (cfg_percdef.trl(:, 1) <= trial_end);
%     
%     percepts(iTrial).trl = cfg_percdef.trl(currPercepts, :);
%     
%     percepts(iTrial).type = dec2bin(cfg_percdef.trl(currPercepts, 4)-1) - '0';
%     
%     percepts(iTrial).start = cfg_percdef.trl(currPercepts, 1);
%     
%     percepts(iTrial).duration = [percepts(iTrial).start(2:end); trial_end] - percepts(iTrial).start;
%     
%     percepts(iTrial).start = (percepts(iTrial).start-trial_start) /all_data.fsample;
%     
%     percepts(iTrial).duration = percepts(iTrial).duration /all_data.fsample;
%     
% end


%% RLS Analysis
rls_data = all_data;
single_rls(1) = all_data;
single_rls(2) = all_data;
msglength = 0;
for iTrial = 1:size(all_data.trial, 2);
    
    fprintf(repmat('\b',1,msglength));
    msg = sprintf('Running RLS - Processing Trial No %d\n', iTrial);
    fprintf(msg);
    msglength = numel(msg);

    colours = {[83 148 255]/255, [255 117 117]/255};
    rls_data.trial{iTrial} = zeros(size(rls_data.trial{iTrial}));

    % Fundamental Stimulation Frequency
    for iFreq = stimulation_frequencies
        j = find(iFreq==stimulation_frequencies);
        cfg_rls.n_cycles = (trial_dur - discard_start) * iFreq;
        cfg_rls.stim_freq = iFreq;
        cfg_rls.channel = focus_electrode;
  
        [single_rls(j).trial{iTrial}, single_rls(j).amp{iTrial}] = rls_slave( cfg_rls, all_data.trial{iTrial} );
        
        rls_data.trial{iTrial} = rls_data.trial{iTrial} + single_rls(j).trial{iTrial};
    end


end

%% Plot the RLS Signal Amplitude Estimate

timeFig = figure;
% scatterFig = figure;
for iTrial = 1:size(all_data.trial, 2);
    
    figure(timeFig);
    subplot(round(trialNum/2), 4, iTrial); hold on;
    
    for iFreq = 1:2;
        
        y{iFreq} = smooth(single_rls(iFreq).amp{iTrial}(focus_electrode, :), 150);        
        plot(y{iFreq}, 'Color', colours{iFreq});
 
    end
    
%     figure(scatterFig);
%     
%     subplot(4, 4, iTrial); hold on; xlim([-4, 4]); ylim([-4, 4]);
%     
%     time_index = single_rls(1).time{1} > 2;
%     
%     xscatt = zscore(single_rls(1).amp{iTrial}(focus_electrode, time_index));
%     yscatt = zscore(single_rls(2).amp{iTrial}(focus_electrode, time_index));
%     
%     scatter(xscatt, yscatt, 1);
end



%% FFT Analysis
% compare the FFT amplitude between RLS estimate and raw data 
cfg_fft = [];
cfg_fft.continuous = 'yes';
cfg_fft.output = 'pow';
cfg_fft.method = 'mtmfft';
cfg_fft.keeptrials = 'yes'; %yes
cfg_fft.foilim = [0 60]; %[0 60]
cfg_fft.tapsmofrq = 0.09;
cfg_fft.channel = 1:32;

freqs_raw = ft_freqanalysis(cfg_fft, all_data);
freqs_rls = ft_freqanalysis(cfg_fft, rls_data);

% figure;
% electrode_to_plot = focus_electrode;
% for iTrial = 1:trialNum 
%     subplot(round(trialNum/4), 5, iTrial); hold on;
%     
%     plot(freqs_raw.freq, squeeze( freqs_raw.powspctrm(iTrial, electrode_to_plot, :) ), 'Color', colours{1},'LineWidth',2);
%     
%     plot(freqs_rls.freq, squeeze( freqs_rls.powspctrm(iTrial, electrode_to_plot, :) ), '--r');
%     xlim([min(stimulation_frequencies)-5, max(stimulation_frequencies)+5]);
% %     ylim([0, 1]);
% %     xlim([0, max([10, stimulation_frequencies])+10]);
% end

% hhh
% return;


%% SNR

for iFreq = 1:length(stimulation_frequencies)
    
    snr(iFreq, 1:32) = ...
        freqs_raw.powspctrm(iTrial, :, find(freqs_raw.freq>stimulation_frequencies(iFreq), 1))...
        ./median(freqs_raw.powspctrm(iTrial, :, (freqs_raw.freq>stimulation_frequencies(iFreq)-3 & freqs_raw.freq<stimulation_frequencies(iFreq)+3)), 3);
    
end


% return;



%% TFR Analysis
% Trial Definition
cfg_tfr = [];
cfg_preproc = [];
discard_start = 0;
[cfg_tfr, cfg_preproc] = configParams(fileID,discard_start,trial_dur,trialNum,electrodes);
all_data = ft_preprocessing(cfg_preproc);

slidingWindowStep = 0.05; %seconds
slidingWindowWidth = 1.2; %window specific for each frequency
maxFreq = 20; %hz %30

cfg_tfr = [];
cfg_tfr.channel = electrodes; %'all' or electrodes
cfg_tfr.foi = 0:(1/30):maxFreq; %frequency of interest %step by 0.25
cfg_tfr.toi = -3:slidingWindowStep:trial_dur; %time of interest
cfg_tfr.keeptrials = 'yes'; % this makes sure all trials are separate

cfg_tfr.method = 'mtmconvol'; % cfg_tfr.method = 'wavelet';
cfg_tfr.taper = 'hanning'; %sliding window shape (tapered)
cfg_tfr.t_ftimwin = slidingWindowWidth*ones(size(cfg_tfr.foi)); % cfg_tfr.width = 16;

tfr_raw = ft_freqanalysis(cfg_tfr, all_data);
 tfr_rls = ft_freqanalysis(cfg_tfr, rls_data);
% 
% for iFreq = 1:2
%  tfr_single_freq{iFreq} = ft_freqanalysis(cfg_tfr, single_rls(iFreq));
% end
cmap = cbrewer('seq', 'PuBu', 200);

yLims = [1 0.8; 1 1; 1 1]; %trial 3
zLims = [0.01 .1; 0.01 .07; 0.01 0.07]; %trial 5
for iTrial = 1:trialNum
    figure;
    for iFreq = 1:length(stimulation_frequencies)
        subplot(2, 1, iFreq);
        cfg_plot = [];
        cfg_plot.channel = freqs_raw.label(focus_electrode);
        cfg_plot.channel = 'all';
        cfg_plot.interactive = 'no';
        cfg_plot.trials = iTrial;
        cfg_plot.xlim = [0 trial_dur];
        cfg_plot.ylim = [stimulation_frequencies(iFreq)-yLims(iFreq,1), stimulation_frequencies(iFreq)+yLims(iFreq,2)];
        %cfg_plot.zlim = [zLims(iFreq,1),zLims(iFreq,2)];
        %cfg_plot.baseline = [-3, -1];
        cfg_plot.baselinetype = 'relative';

        gridxy2([], stimulation_frequencies);

        ft_singleplotTFR(cfg_plot, tfr_rls  );
        %colormap(cmap);
    end
end

return;

%% Some Plotting

figure; 
subplot(3, 1, 1); hold on;
x = tfr_raw.time;
y1 = squeeze(tfr_raw.powspctrm(1, 1, :,:));
y2 = squeeze(tfr_raw.powspctrm(1, 2, :,:));
plot(x, y1, 'Color', colours{1});
plot(x, y2, 'Color', colours{2});

% 
% subplot(3, 1, 2); hold on;
% x = tfr_rls.time;
% y1 = squeeze(tfr_rls.powspctrm(1, 1, :));
% y2 = squeeze(tfr_rls.powspctrm(1, 2, :));
% plot(x, y1, 'Color', colours{1});
% plot(x, y2, 'Color', colours{2});
% 
% subplot(3, 1, 3); hold on;
% x = tfr_single_freq{1}.time;
% y1 = squeeze(tfr_single_freq{1}.powspctrm(1, 1, :));
% y2 = squeeze(tfr_single_freq{2}.powspctrm(1, 2, :));
% plot(x, y1, 'Color', colours{1});
% plot(x, y2, 'Color', colours{2});

% figure;
% hold on;
% for iTrial = 1:16
%     subplot(4, 4, iTrial); hold on;
%     includeTime = tfr_single_freq{1}.time > 0.5;
%     x = squeeze(tfr_single_freq{1}.powspctrm(iTrial, 1, 1, includeTime));
%     y = squeeze(tfr_single_freq{2}.powspctrm(iTrial, 1, 2, includeTime));
%     scatter(x, y);
%     lsline;
% end
% 
% figure;
% hold on;
% for iTrial = 1:16
%     subplot(4, 4, iTrial); hold on;
%     
% %     includeTime = single_rls(j).smooth_time{iTrial} > 1;
%     
%     x = squeeze(single_rls(1).amp{iTrial}(29, :));
%     
% %     y = squeeze(single_rls(2).smooth_amp{iTrial}(29, includeTime));
%     
%     y = single_rls_sub(1).amp{iTrial}(29, :);
%     
%     scatter(x, y);
%     lsline;
% end


%% Amplitude during Percepts

cfg_trl_tfr = [];
figure;

for iTrial = 1:16
    
    % Left Trials
    left_index = percepts(iTrial).type(:, 1) & ~( percepts(iTrial).type(:, 2) | percepts(iTrial).type(:, 3) );
    right_index = percepts(iTrial).type(:, 3) & ~( percepts(iTrial).type(:, 2) | percepts(iTrial).type(:, 1) );
    mid_index = ~( left_index | right_index );
    
    % Empty the time point logicals
    left_time_points = false(size(single_rls(1).time{iTrial}));
    right_time_points = false(size(single_rls(1).time{iTrial}));
    mid_time_points = false(size(single_rls(1).time{iTrial}));
    
    % Go through each type of percept and index the time points within them
    for t = find(left_index)';
        left_time_points = single_rls(1).time{iTrial} >= percepts(iTrial).start(t) & ...
                            single_rls(1).time{iTrial} < percepts(iTrial).start(t) + percepts(iTrial).duration(t) | ...
                            left_time_points;
    end
    for t = find(right_index)';
        right_time_points = single_rls(1).time{iTrial} >= percepts(iTrial).start(t) & ...
                            single_rls(1).time{iTrial} < percepts(iTrial).start(t) + percepts(iTrial).duration(t) | ...
                            right_time_points;
    end
    for t = find(mid_index)';
        mid_time_points = single_rls(1).time{iTrial} >= percepts(iTrial).start(t) & ...
                            single_rls(1).time{iTrial} < percepts(iTrial).start(t) + percepts(iTrial).duration(t) | ...
                            mid_time_points;
    end
    
    % On a trial-by-trial basis, was A or B higher or lower during LEFT or
    % RIGHT
    amp_28_left = mean( single_rls(1).amp{iTrial}(29, left_time_points) ) - mean(single_rls(1).amp{iTrial}(29, :));
    amp_28_right = mean( single_rls(1).amp{iTrial}(29, right_time_points) ) - mean(single_rls(1).amp{iTrial}(29, :));
    amp_36_left = mean( single_rls(2).amp{iTrial}(29, left_time_points) ) - mean(single_rls(2).amp{iTrial}(29, :));
    amp_36_right = mean( single_rls(2).amp{iTrial}(29, right_time_points) ) - mean(single_rls(2).amp{iTrial}(29, :));
    
    subplot(4, 4, iTrial); hold on;
    plot([1, 2], [amp_28_left, amp_28_right], 'Color', colours{1});
    plot([1, 2], [amp_36_left, amp_36_right], 'Color', colours{2});
    
    
end

