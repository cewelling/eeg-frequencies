function [ meanSNR ] = runFFT(  EEGfile, runName, parName, iPar, date, paramsFlag)
%runFFT Summary of this function goes here
%   Detailed explanation goes here

% Ensure path is in correct order
path(genpath('/Users/acspiegel/Documents/MATLAB/Toolboxes/FieldTrip'), path);

%load parameters
analysisParams

%% Update the progress bar
if ~isempty(iPar)
    clc; fprintf('Processing FFT: Subject %2.0f of %2.0f\n', iPar, length(parNums));
else
    clc; fprintf('Processing FFT for group analysis');
end

%% Trial definition and pre-processing

[all_data] = ftConfigParams(parName, runName, EEGfile, numTrials, trialDur, date, [], analFreqs(1) - 2, analFreqs(2) + 2);
[allFreqs_data] = ftConfigParams(parName, runName, EEGfile, numTrials, trialDur, date, [], 2, 70);

%% FFT for each electrode set

analElectrodes = electrodeSet.nums;

cfg_fft = [];
cfg_fft.output = 'pow'; % returns the power-spectra
cfg_fft.method = 'mtmfft'; % analyzes entire spectrum
cfg_fft.foilim = [5, 30]; %[5, 30];
cfg_fft.taper = 'hanning';
% cfg_fft.tapsmofrq = .05; % multitaper spectral smoothing in Hz
cfg_fft.channel = analElectrodes; % electrodes or occipitals
cfg_fft.keeptrials = 'yes'; % return individual trials, not average
% cfg_fft.pad = ceil(max(concatData.time{1,1})); % might be useful
% if we want to compare spectra between groups

freqs = ft_freqanalysis(cfg_fft, all_data); %all_data, concatData, high/lowFreqReportData
% all_data is pre-processed data
allFreqs = ft_freqanalysis(cfg_fft, allFreqs_data);  % FFT with all frequencies

FFTsToSave = [];

for iTrial = 0:size(freqs.powspctrm,1)
    
    if iTrial == 0 % mean across trials
        powspctrm = squeeze(nanmean(freqs.powspctrm,1)); % power per frequency bin for each channel
    else % individual trials
        powspctrm = squeeze(freqs.powspctrm(iTrial,:,:)); % in this trial, power per frequency bin for each channel
        powspctrm_AF = squeeze(allFreqs.powspctrm(iTrial,:,:)); % includes all frequencies
    end
    
    %% SNR Calculation
    
    % calculate SNR for each frequency of interest
    meanSNRs = struct([]);
    for stimFreq = analFreqs;
        % find the index of the FOI
        [~, stimFreqIndex] = min(abs(freqs.freq - stimFreq));
        
        for noiseWindowHalf = noiseWindowHalves;
            
            % identify noise window, exclude the FOI
            noiseIndices1 = [stimFreqIndex - noiseWindowHalf:stimFreqIndex - 1];
            noiseIndices2 = [stimFreqIndex + 1 : stimFreqIndex + noiseWindowHalf];
            noiseIndices = [noiseIndices1 noiseIndices2];
            
            % if only one electrode is used, power is presented in a column
            % instead of a row
            if size(powspctrm, 2) == 1
                powspctrm = powspctrm';
            end
            
            % calculate SNR
            signal = powspctrm(:, stimFreqIndex);
            noise = nanmean(powspctrm(:,noiseIndices), 2);
            SNR = signal./noise;
            
            % Store SNRs for each electrode...note that this only works if
            % we've only chosen one noiseWindow
            if(stimFreq == analFreqs(1))
                SNRlow = SNR;
            else
                SNRhigh = SNR;
            end
            
            newMeanSNR.value = nanmean(SNR); % mean across electrodes
            newMeanSNR.freq = stimFreq;
            newMeanSNR.noiseWindow = noiseWindowHalf*2;
            meanSNRs = [meanSNRs, newMeanSNR];
        end
    end
    SNRs = mean([SNRlow SNRhigh], 2); % mean across frequencies
    
    % Select the maximally responding electrodes
    %          %sort by SNR to find best electrodes
    
    SNRs(isnan(SNRs)) = -inf; % ensure that nans go to end of sorted SNRs
    [sortedSNRs,sortingIndices] = sort(SNRs,'descend');
    
    % order electrode labels from highest to lowest SNR
    maxSNRelecs = [];
    for i = 1:length(analElectrodes)
        maxSNRelecs = [maxSNRelecs analElectrodes(sortingIndices(i))];
    end
    
    % Break down sorted mean SNRS into high freq and low freq SNRs
    lowSNRs = SNRlow(sortingIndices);
    highSNRs = SNRhigh(sortingIndices);
    
    maxSNRs = [];
    lFreqSNRs = []; % SNRs of low frequency band
    hFreqSNRs = []; % SNRs of high frequency band
    
    % electrode indices in first row
    maxSNRs(1,:) = maxSNRelecs;
    lFreqSNRs(1,:) = maxSNRelecs;
    hFreqSNRs(1,:) = maxSNRelecs;
    
    % SNR value in second row
    maxSNRs(2,:) = sortedSNRs;
    lFreqSNRs(2,:) = lowSNRs;
    hFreqSNRs(2,:) = highSNRs;
    
    % Account for mistakes
    if ~isempty(strfind(EEGfile, 'cumulus05_marzRival1'))
        % account for loss of first two trials of this run
        thisTrial = iTrial + 2;
    elseif ~isempty(strfind(EEGfile, 'stratus135_dartSim1'))
        % account for loss of first trial of this run
        thisTrial = iTrial + 1;
    elseif ~isempty(strfind(EEGfile, 'cumulus01_dartSim3'))
        % account for loss of first trial of this run
        thisTrial = iTrial +1;
    elseif ~isempty(strfind(EEGfile, 'cumulus01_marzRival3'))
        % account for loss of first two trials of this run
        thisTrial = iTrial +2;
    else
        thisTrial = iTrial;
    end
    
    if sum(analFreqs == stimFreqs) == 2 % both frequencies are correct
        if iTrial == 0
            save(['pre-processing/highSNRelecs/' parName '_' runName '_'  electrodeSet.name], 'maxSNRs');
        else
            save(['pre-processing/highSNRelecs/' parName '_' runName '_trial' num2str(thisTrial) '_'  electrodeSet.name], 'maxSNRs', 'lFreqSNRs', 'hFreqSNRs');
        end
    elseif sum(analFreqs == imFreqs) == 2
        if iTrial == 0
            save(['pre-processing/highSNRelecs/' parName '_' runName '_'  electrodeSet.name '_IMs'], 'maxSNRs');
        else
            save(['pre-processing/highSNRelecs/' parName '_' runName '_trial' num2str(thisTrial) '_'  electrodeSet.name '_IMs'], 'maxSNRs', 'lFreqSNRs', 'hFreqSNRs');
        end
    end

    % load noisy electrodes list
    load(['pre-processing/badElecs/' parName '_' runName '.mat']);
    %     if isempty(badElecs)
    %         badNums = [];
    %     else
    %         badNums = extractfield(badElecs, 'num');
    %     end
    
    %% Store power spectrum (if SNR is high enough)
  
    if iTrial ~= 0
        
        % load SNRs for principle frequencies
        load(['pre-processing/highSNRelecs/' parName '_' runName '_trial' num2str(thisTrial) '_'  electrodeSet.name '.mat']);
        
        if nanmean(lFreqSNRs(2,lFreqSNRs(1,:)==snrElecs)) >= minSNR && ...
                nanmean(hFreqSNRs(2,hFreqSNRs(1,:)==snrElecs)) >= minSNR % if SNR is high enough
            FFTsToSave = cat(3, FFTsToSave, powspctrm_AF);
        end
    end
    
    %% Plot power spectrum
    if strcmp(FFTplotOrNot, 'yes')
        figure;
        %plot(freqs.freq, powspctrm(analElectrodes,:));
        %plot(freqs.freq, powspctrm(analElectrodes(~ismember(analElectrodes, badElecs)),:)); % exclude noisy electrodes
        %plot(freqs.freq, powspctrm(maxSNRs(1,:),:)); % exclude noisy electrodes
        plot(freqs.freq, powspctrm); % exclude noisy electrodes
        
        vline(analFreqs(1))
        vline(analFreqs(2))
        xlabel('Frequency (Hz)')
        ylabel('Power')
        %plot(freqs.freq, powspctrm(maxIndices(~ismember(maxIndices, badNums)),:)); % exclude noisy electrodes
        %             axis([5 30 0 0.07]) %0.18
        %     gridxy([5, 12]);
        %     xlim([3, 10]);
        %     ylim([0, 15]);
        %     gridxy(stimFreqs);
        box off
        if sum(analFreqs == stimFreqs) == 2 % both frequencies are correct
            if iTrial == 0
                title(['FFT Results: ' parName ', ' runName ', ' electrodeSet.name ])
                figName = [fftSpectPlotDir parName '_' runName '_FFT_' electrodeSet.name '.jpg'];
            else
                title(['FFT Results: ' parName ', ' runName ', trial ' num2str(thisTrial) ', ' electrodeSet.name ])
                figName = [fftSpectPlotDir parName '_' runName '_trial' num2str(thisTrial) '_FFT_' electrodeSet.name '.jpg'];
            end
        elseif sum(analFreqs == imFreqs) == 2 % both frequencies are correct
            if iTrial == 0
                title(['FFT Results: ' parName ', ' runName ', ' electrodeSet.name ])
                figName = [fftSpectPlotDir parName '_' runName '_FFT_' electrodeSet.name '_IMs.jpg'];
            else
                title(['FFT Results: ' parName ', ' runName ', trial ' num2str(thisTrial) ', ' electrodeSet.name ])
                figName = [fftSpectPlotDir parName '_' runName '_trial' num2str(thisTrial) '_FFT_' electrodeSet.name '_IMs.jpg'];
            end
        end
        saveas(gcf,figName,'jpg')
        %             for iTrial = 1:size(freqs.powspctrm, 1)
        %                 figure;
        %                 plot(freqs.freq, squeeze(freqs.powspctrm(iTrial, 5,:)))
        %                 box off
        %                 title(['FFT Results: ' parName ', ' runName ', trial ' num2str(iTrial) ', ' electrodeSets(iElectrodes).name ])
        %                 vline(8.5)
        %                 vline(5.666)
        %             end
        
        % Plot bad electrodes
        %         for iTrial = 1:length(badElecs)
        %             for iElec = 1:length(badElecs{iTrial})
        %                 if ismember(badElecs{iTrial}(iElec).num, analElectrodes)
        %                     figure;
        %                     plot(freqs.freq, squeeze(freqs.powspctrm(iTrial,badElecs{iTrial}(iElec).num,:)));
        %                     box off
        %                     title([parName ', ' runName ', Electrode' num2str(badElecs{iTrial}(iElec).num) ...
        %                         ', Trial: ' num2str(iTrial) ', Variance: ' num2str(badElecs{iTrial}(iElec).var) ...
        %                         ', Mean: ' num2str(badElecs{iTrial}(iElec).mean) ', Kurtosis: ' num2str(badElecs{iTrial}(iElec).kur)])
        %                 end
        %             end
        %         end
        
        % Plot test electrodes
        load(['pre-processing/testElecs/' parName '_' runName '.mat']);
        for iTrial = 1:length(testElecs)
            for iElec = 1:length(testElecs{iTrial})
                if ismember(testElecs{iTrial}(iElec).num, analElectrodes)
                    figure;
                    plot(freqs.freq, squeeze(freqs.powspctrm(iTrial,testElecs{iTrial}(iElec).num,:)));
                    box off
                    title(['FFT Results: ' parName ', ' runName ', Electrode' num2str(testElecs{iTrial}(iElec).num) ...
                        ', Variance: ' num2str(testElecs{iTrial}(iElec).var) ', Kurtosis: ' num2str(testElecs{iTrial}(iElec).kur) ...
                        ', Mean: ' num2str(testElecs{iTrial}(iElec).mean)])
                end
            end
        end
        
      
    end
end

%% Save power spectrum

if ~isempty(FFTsToSave)
    freqAxis = allFreqs.freq;
    for iElec = 1:length(electrodeSet.nums);
        elecFFT = (squeeze(FFTsToSave(iElec,:,:)))';
        saveName = ['FFTs/' parName '_' runName '_' num2str(electrodeSet.nums(iElec))];
        save(saveName, 'elecFFT', 'freqAxis');
    end
end

end

