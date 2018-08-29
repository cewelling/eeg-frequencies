function [maxResponseElecs, maxIndices] = runFFT(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, trialNum, trial_dur,parName,legend,omniData)

%% Analyse Binoc with FFT
% clearvars;
% clc;
%clearvars -except stimulation_frequencies analRun


%% Load default anal params
defaultAnalParams

%%
trialNum = legend.data(analRunNum,3);
trial_dur = legend.data(analRunNum,4);

for groupIndex = 1:1
    for subject = 1:nPars(groupIndex)
        
    %% Update the progress bar
    clc; fprintf('Processing: groupIndex %2.0f of 2, Subject %2.0f of %2.0f\n', groupIndex, subject, nPars(groupIndex));
    fileID = fullfile(file_directory, filenames{groupIndex}(subject).name);


    %% Trial Definition
    cfg_trldef = [];
    cfg_trldef.dataset = fileID;
    cfg_trldef.trialdef.eventtype = 'STATUS';
    cfg_trldef.trialfun = 'ft_trialfun_general';
    cfg_trldef.trialdef.prestim = -discard_start;
    cfg_trldef.trialdef.poststim = trial_dur; 
    cfg_trldef.trialdef.eventvalue = 201:(200+trialNum);
    
    try
        cfg_trldef = ft_definetrial(cfg_trldef);
    catch define_trial_error
        cfg_trldef.trialdef.eventvalue
        cfg_trldef.dataset
        rethrow(define_trial_error);
    end

    %cfg_trldef.trl = remove_overlaps(cfg_trldef); % this script removes overlapping trials in case triggers got confused

    %% Preprocessing

    cfg_preproc = cfg_trldef;
    cfg_preproc.channel = 'all';
    cfg_preproc.continuous = 'yes';
    cfg_preproc.demean    = 'yes'; % perform baseline correction 
    cfg_preproc.detrend = 'yes'; % remove linear trend from the data
    cfg_preproc.reref = 'yes';
    cfg_preproc.refchannel = 'all'; % common average reference

    cfg_preproc.hpfilter = 'yes'; % highpass filter
    cfg_preproc.hpfreq = 2; % highpass frequency in hertz

    all_data = ft_preprocessing(cfg_preproc);
    
    %cfg        = [];
    %cfg.method = 'channel';
    %ft_rejectvisual(cfg, all_data)
    
    % Find and remove artifacts
%     clean_data = removeArtifacts(all_data,cfg_trldef,'yes');
    
    % Concatenate trials
%     concatData = concatenateTrials(all_data);
%     
%     all_data = concatData;
    
    
    % Create data files corresponding to perceptual reports
%     [lowFreqReportData, highFreqReportData] = epochOfInt(omniData, concatData);
    
        for electrodeSet = 1:2

            if electrodeSet == 1
                analElectrodes = electrodes;
                electrodeNames = 'allElectrodes';
            else
%                 analElectrodes = occipitals;
%                 electrodeNames = 'occipitals';
                analElectrodes = 29;
                electrodeNames = 'Oz';
            end

            %% FFT: All electrodes
%             cfg_fft = [];
%             cfg_fft.continuous = 'yes';
%             cfg_fft.output = 'pow';
%             cfg_fft.method = 'mtmfft';
%             cfg_fft.foilim = [5, 30]; %[5, 30];
%             cfg_fft.tapsmofrq = 0.04; %0.09 (about the limit for 12s trials), able to go as low as ~.035
%             cfg_fft.channel = analElectrodes; % electrodes or occipitals
%             cfg_fft.keeptrials = 'no';

            cfg_fft = [];
%             cfg_fft.continuous = 'yes';
            cfg_fft.output = 'pow'; % returns the power-spectra
            cfg_fft.method = 'mtmfft'; % analyzes entire spectrum
            cfg_fft.foilim = [5, 30]; %[5, 30];
            cfg_fft.taper = 'hanning';
%             cfg_fft.tapsmofrq = .05; %0.09 (about the limit for 12s trials), able to go as low as ~.035
            cfg_fft.channel = analElectrodes; % electrodes or occipitals
            cfg_fft.keeptrials = 'no'; % return average, not individual trials
%             cfg_fft.pad = ceil(max(concatData.time{1,1}));
            
            freqs = ft_freqanalysis(cfg_fft, all_data); %all_data, concatData, high/lowFreqReportData
            % all_data is pre-processed data
            groupIndex_freqs{groupIndex, subject} = freqs;
            
            %calculate snr
            dataSNR = squeeze( freqs.powspctrm(:, :) ); %electrodes by frequency
            for noiseWindow = 1:1:19
                %noiseWindow = 4; %hz (units are not integers, units are sampling rate of EEG)
                stimFreqIndices = [max(find(freqs.freq < stimulation_frequencies(1))) max(find(freqs.freq < stimulation_frequencies(2)))];
                for freqIndex = 1:size(stimulation_frequencies,2)
                   signal(:,freqIndex) = dataSNR(:,stimFreqIndices(freqIndex));
                   noise(:,freqIndex) = nanmean(dataSNR(:,stimFreqIndices(freqIndex)-noiseWindow:stimFreqIndices(freqIndex)+noiseWindow),2);
                end
                SNR = signal./noise;
                meanSNR = nanmean(SNR,1) %mean across electrodes for each of the 2 frequencies
            end


            %% Plot the spectrum        
            figure;
            plot(freqs.freq, squeeze( freqs.powspctrm(:, :) ));
%             axis([5 30 0 0.07]) %0.18
            %     gridxy([5, 12]);
            %     xlim([3, 10]);
            %     ylim([0, 15]);
            %     gridxy(stimulation_frequencies);
            box off
            title(['FFT Results: ' parName ', run ' num2str(analRunNum) ' , ' electrodeNames ])
            figName = [fftSpectPlotDir parName '_run' num2str(analRunNum) '_FFT_' electrodeNames '.jpg']
            saveas(gcf,figName,'jpg')
           
            %% Plot the spectrum - means        
%             figure;
%             plot(freqs.freq, mean(squeeze( freqs.powspctrm(:, :) )));
% %             axis([5 30 0 0.07]) %0.18
%             %     gridxy([5, 12]);
%             %     xlim([3, 10]);
%             %     ylim([0, 15]);
%             %     gridxy(stimulation_frequencies);
%             box off
%             title(['FFT Results: ' parName ', run ' num2str(analRunNum) ' , ' electrodeNames ])


            %% Select the maximally responding electrodes
            if electrodeSet == 5 %do this when analyzing occipitals (not used currently)
                
                % find the stim freq indices
                freqIndices = [];
                for i = 1:length(stimulation_frequencies)
                    [c index] = min(abs(freqs.freq-stimulation_frequencies(i)));
                    freqIndices = [freqIndices index];
                end
                
                %take the mean of electrode values at the two stim freq
                %indices
                freqResponse = [freqs.powspctrm(:,freqIndices(1)) freqs.powspctrm(:,freqIndices(2))];
                freqResponse = mean(freqResponse,2);
                
                %sort to find best electrodes
                [sortedElectrodes,sortingIndices] = sort(freqResponse,'descend');
                
                maxResponseElecs = [freqs.label(sortingIndices(1)); freqs.label(sortingIndices(2))];
                maxIndices = [];
                for i = 1:length(maxResponseElecs)
%                     labelMatch = strfind(all_data.label,maxResponseElecs(1,1));
                    maxIndices(i) = find(ismember(all_data.label,maxResponseElecs(i)));
                end
                    
                
            end
                
        end

    end % End of participant loop
end % End of group loop

