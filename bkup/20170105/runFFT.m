function [ meanSNRs ] = runFFT(  EEGfile, runName, parName, iPar, date)
%runFFT Summary of this function goes here
%   Detailed explanation goes here

%load parameters
analysisParams

%% Update the progress bar
clc; fprintf('Processing FFT: Subject %2.0f of %2.0f\n', iPar, length(parNums));

%% Trial definition and pre-processing

[all_data] = ftConfigParams(parName, runName, EEGfile, numTrials, trialDur, date);

%% FFT for each electrode set

for iElectrodes = 1:size(electrodeSets, 2)
    analElectrodes = electrodeSets(iElectrodes).nums;
    
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
    
    powspctrm = squeeze(nanmean(freqs.powspctrm,1)); % power per frequency bin for each channel
    
    %% SNR Calculation
    
    % calculate SNR for each frequency of interest
    meanSNRs = struct([]);
    for stimFreq = stimFreqs;
        % find the index of the FOI
        [~, stimFreqIndex] = min(abs(freqs.freq - stimFreq));
        
        for noiseWindowHalf = noiseWindowHalves;
            
            % identify noise window, exclude the FOI
            noiseIndices1 = [stimFreqIndex - noiseWindowHalf:stimFreqIndex - 1];
            noiseIndices2 = [stimFreqIndex + 1 : stimFreqIndex + noiseWindowHalf];
            noiseIndices = [noiseIndices1 noiseIndices2];
            
            signal = powspctrm(:, stimFreqIndex);
            noise = nanmean(powspctrm(:,noiseIndices), 2);
            SNR = signal./noise;
            
            % Store SNRs for each electrode...note that this only works if
            % we've only chosen one noiseWindow
            if(stimFreq < 6)
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
    %save ([indicesDir 'SNR/' 'MarzAllElectrodes_' parName '_run' num2str(currentRun)], 'SNRs');
    save ([indicesDir 'SNR/' parName '_' runName '_' electrodeSets(iElectrodes).name], 'meanSNRs');
    
    % Select the maximally responding electrodes
    %          %sort by SNR to find best electrodes
    [~,sortingIndices] = sort(SNRs,'descend');
    
    maxSNRelecs = [];
    for i = 1:numElecs
        maxSNRelecs = [maxSNRelecs; freqs.label(sortingIndices(i))];
    end
    maxIndices = [];
    for i = 1:length(maxSNRelecs)
        %                     labelMatch = strfind(all_data.label,maxSNRelecs(1,1));
        maxIndices(i) = find(ismember(freqs.label,maxSNRelecs(i)));
    end
    
    save(['highSNRelecs/' parName '_' runName '_'  num2str(numElecs) 'elecs'], 'maxIndices');
    
    % load noisy electrodes list
    load(['badElecs/' parName '_' runName '.mat']);
%     if isempty(badElecs)
%         badNums = [];
%     else
%         badNums = extractfield(badElecs, 'num');
%     end
    
    %% Plot power spectrum
    if strcmp(FFTplotOrNot, 'yes')
        figure;
        %plot(freqs.freq, powspctrm(analElectrodes,:)); 
        %plot(freqs.freq, powspctrm(analElectrodes(~ismember(analElectrodes, badElecs)),:)); % exclude noisy electrodes      
        plot(freqs.freq, powspctrm(maxIndices,:)); % exclude noisy electrodes
        %plot(freqs.freq, powspctrm(maxIndices(~ismember(maxIndices, badNums)),:)); % exclude noisy electrodes
        %             axis([5 30 0 0.07]) %0.18
        %     gridxy([5, 12]);
        %     xlim([3, 10]);
        %     ylim([0, 15]);
        %     gridxy(stimFreqs);
        box off
        title(['FFT Results: ' parName ', ' runName ', ' electrodeSets(iElectrodes).name ])
        figName = [fftSpectPlotDir parName '_' runName '_FFT_' electrodeSets(iElectrodes).name '.jpg'];
        saveas(gcf,figName,'jpg')
        
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
        load(['testElecs/' parName '_' runName '.mat']);
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
        
        %         for i = 1:length(badElecs)
        %             if ismember(badElecs(i).num, analElectrodes)
        %                 figure;
        %                 plot(freqs.freq, powspctrm(intersect(analElectrodes,badElecs(i).num),:));
        %                 %plot(freqs.freq, powspctrm(analElectrodes(~ismember(analElectrodes, badElecs)),:)); % exclude noisy electrodes
        %                 %plot(freqs.freq, powspctrm(maxIndices,:)); % exclude noisy electrodes
        %                 %plot(freqs.freq, powspctrm(maxIndices(~ismember(maxIndices, badElecs)),:)); % exclude noisy electrodes
        %                 %             axis([5 30 0 0.07]) %0.18
        %                 %     gridxy([5, 12]);
        %                 %     xlim([3, 10]);
        %                 %     ylim([0, 15]);
        %                 %     gridxy(stimFreqs);
        %                 box off
        %                 title(['FFT Results: ' parName ', ' runName ', Electrode' num2str(badElecs(i).num) ...
        %                     ', Variance: ' num2str(badElecs(i).var) ', Kurtosis: ' num2str(badElecs(i).kur)])
        %             end
        %         end
    end  
end
end

