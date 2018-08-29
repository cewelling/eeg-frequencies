clearvars;
clc;

%% Prompt for Run/Trial(s) of Interest

prompt = {'Which testing day?','Run Number?'};
dlg_title = 'Run of Interest';
num_lines = 1;
defaultans = {'6','3'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

% Add button to select which plots to view


%% Analyse Binoc with FFT


%% File Names

n{1} = 1;
analRunNum = str2double(answer(2));

%% Analysis Variables
% trial_dur = 30; %JL 12; %JF
discard_start = 0.5; % how much time should be cut off at beginning
occipitals = [20 26 27 28 29 30 31];
electrodes = occipitals; %[1:32]; %[27, 29, 30,28]; %occipitals

%stimulation_frequencies = [28.33]; %, 21.25];
% intermodulation_frequencies = [21.6, 43.2];
% intermodulation_frequencies = 21.6;

% Dates & Day #: day2 - 3_11; day3 - 3_14; day4 - 3_15; day5 - 3_17; day6 -
% 3_21
datesDays = {'0';'3_11';'3_14';'3_15';'3_17';'3_21'};
date = char(datesDays(str2double(answer(1))));

fileNames = dir(['../../runScripts/Results/caroline*' date '*']);
fileNames = {fileNames.name};
filenames{1}(1).name = char(fileNames(analRunNum));
legend = importdata(['pilotlegends/pilot-day' char(answer(1)) '-caroline.txt']);

stimulation_frequencies = legend.data(analRunNum,1:2);
trialNum = legend.data(analRunNum,3);
trial_dur = legend.data(analRunNum,4);


for group = 1:1
    for subject = 1:n{group}
    %% Update the progress bar, load psychometric data
    clc; fprintf('Processing: Group %2.0f of 2, Subject %2.0f of %2.0f\n', group, subject, n{group});
    fileID = fullfile(file_directory, filenames{group}(subject).name);


    %% Trial Definition
    cfg_trldef = [];
    cfg_trldef.dataset = fileID;
    cfg_trldef.trialdef.eventtype = 'STATUS';
    cfg_trldef.trialfun = 'ft_trialfun_general';
    cfg_trldef.trialdef.prestim = -discard_start;
    cfg_trldef.trialdef.poststim = trial_dur; %actual length of trial 12 s
    cfg_trldef.trialdef.eventvalue = 201:(200+trialNum);

    try
        cfg_trldef = ft_definetrial(cfg_trldef);
    catch define_trial_error
        cfg_trldef.trialdef.eventvalue
        cfg_trldef.dataset
        rethrow(define_trial_error);
    end

    %     cfg_trldef.trl = remove_overlaps(cfg_trldef); % this script removes overlapping trials in case triggers got confused

    %% Preprocessing

    cfg_preproc = cfg_trldef;
    cfg_preproc.channel = 'all';
    cfg_preproc.continuous = 'yes';
    cfg_preproc.demean    = 'yes';
    cfg_preproc.detrend = 'yes';
    cfg_preproc.reref = 'yes';
    cfg_preproc.refchannel = 'all';

    cfg_preproc.hpfilter = 'yes';
    cfg_preproc.hpfreq = 2;

    all_data = ft_preprocessing(cfg_preproc);
    
    %cfg        = [];
    %cfg.method = 'channel';
    %ft_rejectvisual(cfg, all_data)


    %% FFT
    cfg_fft = [];
    cfg_fft.continuous = 'yes';
    cfg_fft.output = 'pow';
    cfg_fft.method = 'mtmfft';
    cfg_fft.foilim = [5, 30]; %[5, 30];
    cfg_fft.tapsmofrq = 0.09;
    cfg_fft.channel = electrodes; % electrodes or occipitals
    cfg_fft.keeptrials = 'no';

    freqs = ft_freqanalysis(cfg_fft, all_data);
    group_freqs{group, subject} = freqs;


    %% Plot the spectrum

    figure;
    for iTrial = 4%1:1;

    %     subplot(4, 4, iTrial);

        plot(freqs.freq, squeeze( freqs.powspctrm(:, :) ));
        axis([5 30 0 0.5]) %0.18
    %     gridxy([5, 12]);

    %     xlim([3, 10]);
    %     ylim([0, 15]);

    %     gridxy(stimulation_frequencies);
    end


    end % End of participant loop
end % End of group loop

%% RLS
%% Analysis Variables

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

figure;
electrode_to_plot = focus_electrode;
for iTrial = 1:trialNum
    subplot(round(trialNum/4), 5, iTrial); hold on;
    
    plot(freqs_raw.freq, squeeze( freqs_raw.powspctrm(iTrial, electrode_to_plot, :) ), 'Color', colours{1},'LineWidth',2);
    
    plot(freqs_rls.freq, squeeze( freqs_rls.powspctrm(iTrial, electrode_to_plot, :) ), '--r');
    xlim([min(stimulation_frequencies)-5, max(stimulation_frequencies)+5]);
%     ylim([0, 1]);
%     xlim([0, max([10, stimulation_frequencies])+10]);
end

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
maxFreq = 30; %hz

cfg_tfr = [];
cfg_tfr.channel = electrodes; %'all' or electrodes
cfg_tfr.foi = 0:0.25:maxFreq; %frequency of interest
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

%% Behavrioral Report

datafile = ['21-Mar-2016_rivalData_caroline0' char(answer(2)) '_day' char(answer(1)) '.txt']; %good one: '14-Mar-2016_rivalData_.txt';
dataDir = '../../runScripts/Results/';
data = load([dataDir datafile]);

upArrow = 38;
leftArrow = 37;
rightArrow = 39;

for trialIndex = 1:5 %1:4
    
   trialData = data(find(data(:,1) == trialIndex),:); 
   
   xdata = 1:length(trialData);
   ydata = trialData(:,2);
   timedata = trialData(:,3);
   lastSecond = round(timedata(end));
   timeMultiplier = lastSecond/xdata(end);
   xdata = xdata * timeMultiplier;
   
   figure
   plot(xdata,ydata)
   box off
   axis([0 lastSecond min(ydata) max(ydata)])

   
   numPoints = ones(1,length(xdata));
   colorData = repmat([0 0 0]',1,length(xdata));
   
   figure
   for index = 1:length(xdata)       
        if ydata(index) == upArrow
            scatter(xdata(index),2,'r')
        elseif ydata(index) == leftArrow
            scatter(xdata(index),2,'b')
        elseif ydata(index) == rightArrow
            scatter(xdata(index),2,'g')
        else
            scatter(xdata(index),2,'k')
        end
        hold on
   end
   axis([0 lastSecond 1.5 2.5])
   box off
   %xaxis = timedata
    
    
end



%% parse rivalry data
cOmniData = data;
cOmniData = cOmniData(find(cOmniData(:,2) > 0),:); %remove times with no button press
minEventDur = 0.25; % Have to have pressed the button for 250 ms to count as a perceptual state

omniStop = []; omniStart = []; omniTime = []; omniKey = [];  omniRow = []; omniTrial = []; omniFreq = [];
for trialIndex = 1:max(cOmniData(:,1))

    cData = cOmniData(find(cOmniData(:,1) == trialIndex),:);
    lastSecond = round(cData(end,3));
    
    %make a time vector that is continuous throughout rivalry experiment
    cData(:,6) = cData(:,3) + (trialIndex-1)*lastSecond; 
    
    %start trial on first state press
    cData = cData(min([min(find(cData(:,2)==leftArrow))],[min(find(cData(:,2)==rightArrow))]):end,:);

    elapsedRows = 0;
    for rowIndex = 1:length(cData)-1

        if cData(rowIndex+1,2) == cData(rowIndex,2)
            elapsedRows = elapsedRows + 1;
        else

            keyIndex = cData(rowIndex,2);
            stopTime = cData(rowIndex,3);
            startTime = cData(rowIndex - elapsedRows,3);
            elapsedTime = stopTime - startTime;
            
            %was the person reporting the fast or slow frequency? (
            freqIndex = 0; %mixed percepts
            if cData(rowIndex,4) == 1; %fast on left
                if keyIndex == leftArrow
                    freqIndex = 1; %fast frequency
                elseif keyIndex == rightArrow
                    freqIndex = -1; %slow frequency
                end
            elseif cData(rowIndex,4) == 2; %fast on right
                if keyIndex == leftArrow
                    freqIndex = -1; %slow frequency
                elseif keyIndex == rightArrow
                    freqIndex = 1; %fast frequency
                end
            end
            
           
            omniFreq = [omniFreq; freqIndex]; %which frequency was seen
            omniStop = [omniStop; stopTime];  %start time of perceptual episode
            omniStart = [omniStart; startTime]; %end time of perceptual episode
            omniTime = [omniTime; elapsedTime]; %duration of perceptual episode
            omniKey = [omniKey; keyIndex]; %which key was pressed
            omniTrial = [omniTrial; trialIndex]; %which rivalry trial we're on
            omniRow = [omniRow; elapsedRows]; 

            elapsedRows = 0;

            clear stop* start* elapsedTime key*                
        end
    end

end

                    %final data
omniData = [omniTrial, omniStart, omniStop, omniTime, omniRow/3, omniKey, omniFreq];
omniData = omniData(find(omniData(:,4) > minEventDur),:)

%% Frequency Power During Behavioral Report

% trialOfInt = 4;

freq1 = 17;
freq2 = 21.25;

pressLatency = 0.4;
% pressLatency = 0.8;
%  pressLatency = 0;

%% format freq data

omniFreqs = [];
omniDiffs = [];
for trialOfInt = 1:5

    spectrum = tfr_rls.powspctrm;
    trialSpect = squeeze(spectrum(trialOfInt,:,:,:));
    meanTrialSpect = squeeze(mean(trialSpect,1));

    timeInd = tfr_rls.time;
    freqInd = tfr_rls.freq;

    %broadband
    % freq1channels = meanTrialSpect(67:69,:); %17 Hz
    % freq2channels = meanTrialSpect(84:86,:); %21 Hz

    freq1channels = meanTrialSpect(56:58,:); %14 Hz
    freq2channels = meanTrialSpect(67:69,:); %17 Hz

    %normalize: should this be done to the max, median, or mean?
    for i = 1:size(freq1channels,1) 
        freq1channels(i,:)  = freq1channels(i,:)/nanmax(freq1channels(i,:));
    end
    for i = 1:size(freq2channels,1) 
        freq2channels(i,:)  = freq2channels(i,:)/nanmax(freq2channels(i,:));
    end


    %narrowband
    % freq1channels = meanTrialSpect(68,:);
    % freq2channels = meanTrialSpect(85,:);

    buttonPresses = omniData;

    trialPresses = [];
    for i = 1:size(buttonPresses,1)
        if buttonPresses(i,1) == trialOfInt
            trialPresses = [trialPresses; buttonPresses(i,:)];
        end
    end

    %convert epoch start-end into time indices
    lowCount = 1;
    for i = 1:size(trialPresses,1)
        if trialPresses(i,7) == -1
            thisEpoch = find(trialPresses(i,2) < timeInd & timeInd < trialPresses(i,3));
            thisEpoch = find((trialPresses(i,2)-pressLatency) < timeInd & timeInd < (trialPresses(i,3)-pressLatency));
            lowIndices(lowCount,1:2) = [thisEpoch(1) thisEpoch(end)];
            lowCount = lowCount + 1;
        end
    end

    % lowIndices = find(23 < timeInd & timeInd < 24); %4-5
    highCount = 1;
    for i = 1:size(trialPresses,1)
        if trialPresses(i,7) == 1
            thisEpoch = find(trialPresses(i,2) < timeInd & timeInd < trialPresses(i,3));
            thisEpoch = find((trialPresses(i,2)-pressLatency) < timeInd & timeInd < (trialPresses(i,3)-pressLatency));
            highIndices(highCount,1:2) = [thisEpoch(1) thisEpoch(end)];
            highCount = highCount + 1;
        end
    end


    % highIndices = find(20 < timeInd & timeInd < 20.5); %14-15

    %freq channel vals during low/high

    freq1low = nanmean(nanmean(freq1channels(:,lowIndices)));
    freq1high = nanmean(nanmean(freq1channels(:,highIndices)));

    freq2low = nanmean(nanmean(freq2channels(:,lowIndices)));
    freq2high = nanmean(nanmean(freq2channels(:,highIndices)));

%     % plot it
%     figure;
%     bar([freq1low freq1high freq2low freq2high])
%     title(['Averaged across events, trial' num2str(trialOfInt)])


    %% 
    for i = 1:length(lowIndices)
        meanFreq1ChannelsLow(i,:) = nanmean(freq1channels(:,lowIndices(i,1):lowIndices(i,2)),2);
        meanFreq2ChannelsLow(i,:) = nanmean(freq2channels(:,lowIndices(i,1):lowIndices(i,2)),2);
        meanDiffLow(i,:) =  meanFreq1ChannelsLow(i,:) - meanFreq2ChannelsLow(i,:);
    end

    for i = 1:length(highIndices)
        meanFreq1ChannelsHigh(i,:) = nanmean(freq1channels(:,highIndices(i,1):highIndices(i,2)),2);
        meanFreq2ChannelsHigh(i,:) = nanmean(freq2channels(:,highIndices(i,1):highIndices(i,2)),2);
        meanDiffHigh(i,:) =  meanFreq2ChannelsHigh(i,:) - meanFreq1ChannelsHigh(i,:);
    end


    freq1low = nanmean(nanmean(meanFreq1ChannelsLow));
    freq1high = mean(mean(meanFreq1ChannelsHigh));

    freq2low = nanmean(nanmean(meanFreq2ChannelsLow));
    freq2high = nanmean(nanmean(meanFreq2ChannelsHigh));
    
    diff1 = nanmean(nanmean(meanDiffHigh));
    diff2 = nanmean(nanmean(meanDiffLow));

    % plot it
    figure;
    bar([freq1low freq1high freq2low freq2high])
    title(['Averaged across events, trial' num2str(trialOfInt)])

    % plot it
    figure;
    bar([diff1 diff2])
    title(['Difference Score: trial' num2str(trialOfInt)])

    omniFreqs(trialOfInt,:) = [freq1low freq1high freq2low freq2high];
    omniDiffs(trialOfInt,:) = [diff1 diff2];

end

%Plot Average Data: Bars
runSum  = mean(omniFreqs);
stds = std(omniFreqs);
figure
hold on
bar(1:4,runSum)
errorbar(1:4,runSum,stds,'.')
title(['Across Trials'])

%Plot Average Data: Difference Scores
runSum  = mean(omniDiffs);
stds = std(omniDiffs);
figure
hold on
bar(1:2,runSum)
errorbar(1:2,runSum,stds,'.')
title(['Across Trials'])
