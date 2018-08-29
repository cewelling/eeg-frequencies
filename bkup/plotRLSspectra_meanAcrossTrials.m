function plotRLSspectra_meanAcrossTrials(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, trialNum, trial_dur,parName,legend,omniData)


%% Analysis Variables
for electrodeSet = 2 %1:2
    
    %% Load default anal params
    defaultAnalParams
    
    if electrodeSet == 1
        analElectrodes = electrodes;
        electrodeNames = 'allElectrodes';
    else
        analElectrodes = occipitals;
        electrodeNames = 'occipitals';
    end


    %% Practice Analysis with one subject
    fileID = fullfile(file_directory, filenames{1}(1).name);

    %% load config params
    [cfg_trldef, cfg_preproc] = configParams(fileID,discard_start,trial_dur,trialNum,analElectrodes)

    %% preprocess data
    all_data = ft_preprocessing(cfg_preproc);
    
    %% Set up nonConcat
    
    % Create data files corresponding to perceptual reports, keeping trials
    % seperate
    [lowFreqReportData, highFreqReportData] = epochOfInt('no',omniData, all_data, pressLatency,'no',0.4);

    
    loopCount = 1;
    for trial = 1:length(all_data.trial) % for each trial in a run session
        for focus_electrode = analElectrodes % for each electrode being considered
            for dataSet = 1:3 %loop through each of the concatenated data sets (complete, low frequency report only, high frequency report only), results in a 3-plot subplot

            %select the data structure to work with
            if dataSet == 1
                dataOfInt = all_data;
                electrodeNames = 'all time points';
            elseif dataSet == 2
                dataOfInt = lowFreqReportData;
                electrodeNames = 'low frequency report time points';
            elseif dataSet == 3
                dataOfInt = highFreqReportData;
                electrodeNames = 'high frequency report time points';
            end

            trial_dur = ceil(max(dataOfInt.time{1,trial}));

            %create structure for this trial
            trial_data = [];
            trial_data.label = dataOfInt.label;
            trial_data.time{1,1} = dataOfInt.time{1,trial};
            trial_data.trial{1,1} = dataOfInt.trial{1,trial};
            trial_data.fsample = dataOfInt.fsample;

            %% RLS Analysis
            rls_data = trial_data;
            single_rls(1) = trial_data;
            single_rls(2) = trial_data;
            msglength = 0;
            for iTrial = 1:size(trial_data.trial, 2); %changed 1 to 2:

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

                    [single_rls(j).trial{iTrial}, single_rls(j).amp{iTrial}] = rls_slave( cfg_rls, trial_data.trial{iTrial} );

                    rls_data.trial{iTrial} = rls_data.trial{iTrial} + single_rls(j).trial{iTrial};
                end


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
            cfg_fft.channel = 1:32; % 1:32

            freqs_raw = ft_freqanalysis(cfg_fft, trial_data);
            freqs_rls = ft_freqanalysis(cfg_fft, rls_data);

%             if dataSet == 1
%                 figure;
%             end
            
%             Plotting each electrode - skip this
%             electrode_to_plot = focus_electrode;
%             for iTrial = 1
%                 subplot(2, 3, dataSet+3);
%                 plot(freqs_raw.freq, squeeze( freqs_raw.powspctrm(iTrial, electrode_to_plot, :) ));
%                 xlim([min(stimulation_frequencies)-5, max(stimulation_frequencies)+5]);
%                 
%                 subplot( 2, 3, dataSet); hold on;
%                 plot(freqs_rls.freq, squeeze( freqs_rls.powspctrm(iTrial, electrode_to_plot, :) ));
%                 xlim([min(stimulation_frequencies)-5, max(stimulation_frequencies)+5]);
%             end
% 
%             box off
%             if dataSet == 2
%                 title(['RLS Fit Results: ' parName ', run ' num2str(analRunNum) ' , electrode number ' num2str(electrode_to_plot)])
%             end
%             figName = [plotDir parName '_run' num2str(analRunNum) '_RLSfits_' electrodeNames '.jpg']
        %     saveas(gcf,figName,'jpg')

            if dataSet == 1
                allTimeRLS(loopCount,:) = squeeze( freqs_rls.powspctrm(iTrial, focus_electrode, :));
                allTimeRaw(loopCount,:) = squeeze( freqs_raw.powspctrm(iTrial, focus_electrode, :));
            elseif dataSet == 2
                lowTimeRLS(loopCount,:) = squeeze( freqs_rls.powspctrm(iTrial, focus_electrode, :));
                lowTimeRaw(loopCount,:) = squeeze( freqs_raw.powspctrm(iTrial, focus_electrode, :));
            elseif dataSet == 3
                highTimeRLS(loopCount,:) = squeeze( freqs_rls.powspctrm(iTrial, focus_electrode, :));
                highTimeRaw(loopCount,:) = squeeze( freqs_raw.powspctrm(iTrial, focus_electrode, :));
            end

            clear rls_data single_rls stimFreq

            end
            loopCount = loopCount + 1;
        end

        allTimeTrialRLS(trial,:) = mean(allTimeRLS);
        lowTimeTrialRLS(trial,:) = mean(lowTimeRLS);
        highTimeTrialRLS(trial,:) = mean(highTimeRLS);
        
        allTimeTrialRaw(trial,:) = mean(allTimeRaw);
        lowTimeTrialRaw(trial,:) = mean(lowTimeRaw);
        highTimeTrialRaw(trial,:) = mean(highTimeRaw);

        loopCount = 1;

        clear concatOccips lowFreqOccips highFreqOccips

    end
    
    % Determine RLS Peak Differences
    
    lowDifferences = [];
    highDifferences = [];
    
    for trial = 1:size(allTimeTrialRLS,1)
        lowTrialDiff = max(lowTimeTrialRLS(trial,420:430)) - max(highTimeTrialRLS(trial,420:430));
        lowDifferences = [lowDifferences; lowTrialDiff];
        highTrialDiff = max(highTimeTrialRLS(trial,505:515)) - max(lowTimeTrialRLS(trial,505:515));
        highDifferences = [highDifferences; highTrialDiff];
    end
    
    lowMinusHigh = mean(lowDifferences);
    lowMinusHighError = ste(lowDifferences);
    
    highMinusLow = mean(highDifferences);
    highMinusLowError = ste(highDifferences);

    % Plot the means with error

    %RLS
    errorallTimeTrialRLS = ste(allTimeTrialRLS);
    meanallTimeTrialRLS = mean(allTimeTrialRLS);

    errorlowTimeTrialRLS = ste(lowTimeTrialRLS);
    meanlowTimeTrialRLS = mean(lowTimeTrialRLS);

    errorhighTimeTrialRLS = ste(highTimeTrialRLS);
    meanhighTimeTrialRLS = mean(highTimeTrialRLS);
    
    %Raw
    errorallTimeTrialRaw = ste(allTimeTrialRaw);
    meanallTimeTrialRaw = mean(allTimeTrialRaw);

    errorlowTimeTrialRaw = ste(lowTimeTrialRaw);
    meanlowTimeTrialRaw = mean(lowTimeTrialRaw);

    errorhighTimeTrialRaw = ste(highTimeTrialRaw);
    meanhighTimeTrialRaw = mean(highTimeTrialRaw);
    

    figure
    subplot(2,3,1);
    mseb(freqs_rls.freq,meanallTimeTrialRLS,errorallTimeTrialRLS)
    xlim([min(stimulation_frequencies)-2, max(stimulation_frequencies)+2]);
    title(['RLS Fit Results: ' parName ', run ' num2str(analRunNum) ' , all time points'])

    subplot(2,3,2);
    mseb(freqs_rls.freq,meanlowTimeTrialRLS,errorlowTimeTrialRLS)
    xlim([min(stimulation_frequencies)-2, max(stimulation_frequencies)+2]);
    title(['RLS Fit Results: ' parName ', run ' num2str(analRunNum) ' , low frequency report time points'])

    subplot(2,3,3);
    mseb(freqs_rls.freq,meanhighTimeTrialRLS,errorhighTimeTrialRLS)
    xlim([min(stimulation_frequencies)-2, max(stimulation_frequencies)+2]);
    title(['RLS Fit Results: ' parName ', run ' num2str(analRunNum) ' , high frequency report time points'])
    
    subplot(2,3,4);
    mseb(freqs_raw.freq,meanallTimeTrialRaw,errorallTimeTrialRaw)
    xlim([min(stimulation_frequencies)-2, max(stimulation_frequencies)+2]);
%     title(['RLS Fit Results: ' parName ', run ' num2str(analRunNum) ' , all time points'])

    subplot(2,3,5);
    mseb(freqs_raw.freq,meanlowTimeTrialRaw,errorlowTimeTrialRaw)
    xlim([min(stimulation_frequencies)-2, max(stimulation_frequencies)+2]);
%     title(['RLS Fit Results: ' parName ', run ' num2str(analRunNum) ' , low frequency report time points'])

    subplot(2,3,6);
    mseb(freqs_raw.freq,meanhighTimeTrialRaw,errorhighTimeTrialRaw)
    xlim([min(stimulation_frequencies)-2, max(stimulation_frequencies)+2]);
%     title(['RLS Fit Results: ' parName ', run ' num2str(analRunNum) ' , high frequency report time points'])

    figure
    subplot(1,2,1)
    title(['Low vs. High Report: RLS, ' parName ', run ' num2str(analRunNum) ', delay: ' num2str(pressLatency)]) 
    mseb(freqs_rls.freq,[meanlowTimeTrialRLS; meanhighTimeTrialRLS],[errorlowTimeTrialRLS; errorhighTimeTrialRLS])
    xlim([min(stimulation_frequencies)-2, max(stimulation_frequencies)+2]);
%     legend('14 Hz Report','17 Hz Report');

    subplot(1,2,2)
    hold on
    title('Differences')
    bar(1:2,[lowMinusHigh highMinusLow])
    Labels = {'14.166 Hz', '17 Hz'};
    set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
    errorbar(1:2,[lowMinusHigh highMinusLow],[lowMinusHighError highMinusLowError],'.')

%     figure
% %     title(['Low vs. High Report: RLS, delay:' num2str(pressLatency
%     mseb(freqs_raw.freq,[meanlowTimeTrialRaw; meanhighTimeTrialRaw],[errorlowTimeTrialRaw; errorhighTimeTrialRaw])
%     xlim([min(stimulation_frequencies)-2, max(stimulation_frequencies)+2]);

  clearvars -except *frequencies  analRunNum  filenames  trialNum  trial_dur parName legend
end
