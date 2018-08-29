function plotRLSspectra(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, trialNum, trial_dur,parName,legend,omniData)

%% Clear everything and establish where data is
% clearvars;
% clc;


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

%% Set up

    % Find and remove artifacts
    clean_data = removeArtifacts(cfg_trldef,'no');

    % Concatenate trials
    concatData = concatenateTrials(all_data);
    
    % Create data files corresponding to perceptual reports
    for latency =  pressLatency; % 0:0.2:1.6; loop through possible button lags to look for effect
    [lowFreqReportData, highFreqReportData] = epochOfInt('yes',omniData, concatData, latency);
    
    loopCount = 1;
    for focus_electrode = analElectrodes %occipitals %[9 10 11 13 14 15 16 25] 'effect region'
        for dataSet = 1:3 %loop through each of the concatenated data sets (complete, low frequency report only, high frequency report only), results in a 3-plot subplot

            if dataSet == 1
                all_data = concatData;
                electrodeNames = 'all time points';
%                 stimFreq = stimulation_frequencies; %stimFreq select to
%                 fit the RLS model to only the frequency reported
            elseif dataSet == 2
                all_data = lowFreqReportData;
                electrodeNames = 'low frequency report time points';
%                 stimFreq = min(stimulation_frequencies);
            elseif dataSet == 3
                all_data = highFreqReportData;
                electrodeNames = 'high frequency report time pooints';
%                 stimFreq = max(stimulation_frequencies);
            end

            %create an additional channel, that is the mean across the occipital
            %electrodes of interest

            trial_dur = ceil(max(all_data.time{1,1}));
            trialNum = 1;

        %     averagedOccipitals = []; % use this to get standard deviation
        %     for i = 1:length(occipitals)
        %         averagedOccipitals(i,:) = all_data.trial{1,1}(occipitals(i),:);
        %     end
        %     
        %     all_data.trial{1,1}(33,:) = mean(averagedOccipitals,1);
        %     all_data.label{33,1} = 'Occ';

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
                for iFreq = stimulation_frequencies %= stimFreq % set to stimFreq in order to model for only the frequency being reported
                    j = find(iFreq==stimulation_frequencies);
                    cfg_rls.n_cycles = (trial_dur - discard_start) * iFreq;
                    cfg_rls.stim_freq = iFreq;
                    cfg_rls.channel = focus_electrode; %focus_electrode (29)

                    [single_rls(j).trial{iTrial}, single_rls(j).amp{iTrial}] = rls_slave( cfg_rls, all_data.trial{iTrial} );

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

            freqs_raw = ft_freqanalysis(cfg_fft, all_data);
            freqs_rls = ft_freqanalysis(cfg_fft, rls_data);

            if dataSet == 1
                figure;
            end

            electrode_to_plot = focus_electrode; %focus_electrode
            for iTrial = 1:trialNum
                subplot( 1, 3, dataSet); hold on; % changed round to ceil: round(trialNum/4), 5

%                 plot(freqs_raw.freq, squeeze( freqs_raw.powspctrm(iTrial, electrode_to_plot, :) ), 'Color', colours{1},'LineWidth',2);

%                 plot(freqs_rls.freq, squeeze( freqs_rls.powspctrm(iTrial, electrode_to_plot, :) ), '--r');
                plot(freqs_rls.freq, squeeze( freqs_rls.powspctrm(iTrial, electrode_to_plot, :) ));
                xlim([min(stimulation_frequencies)-5, max(stimulation_frequencies)+5]);
            %     ylim([0, 1]);
            %     xlim([0, max([10, stimulation_frequencies])+10]);
            end

            box off
            if dataSet == 2
                title(['RLS Fit Results: ' parName ', run ' num2str(analRunNum) ' , electrode number ' num2str(electrode_to_plot)])
            end
            figName = [plotDir parName '_run' num2str(analRunNum) '_RLSfits_' electrodeNames '.jpg']
        %     saveas(gcf,figName,'jpg')
        
            if dataSet == 1
                concatOccips(loopCount,:) = squeeze( freqs_rls.powspctrm(iTrial, electrode_to_plot, :));
            elseif dataSet == 2
                lowFreqOccips(loopCount,:) = squeeze( freqs_rls.powspctrm(iTrial, electrode_to_plot, :));
            elseif dataSet == 3
                highFreqOccips(loopCount,:) = squeeze( freqs_rls.powspctrm(iTrial, electrode_to_plot, :));
            end

            clear rls_data single_rls stimFreq

        end
        loopCount = loopCount + 1;
    end
    
    % Plot the means with error
    
    errorConcatOccips = ste(concatOccips);
    meanConcatOccips = mean(concatOccips);
    
    errorlowFreqOccips = ste(lowFreqOccips);
    meanlowFreqOccips = mean(lowFreqOccips);
    
    errorhighFreqOccips = ste(highFreqOccips);
    meanhighFreqOccips = mean(highFreqOccips);
    
    figure
    subplot(1,3,1);
    mseb(freqs_rls.freq,meanConcatOccips,errorConcatOccips)
    xlim([min(stimulation_frequencies)-2, max(stimulation_frequencies)+2]);
    title(['RLS Fit Results: ' parName ', run ' num2str(analRunNum) ' , all time points'])
    
    subplot(1,3,2);
    mseb(freqs_rls.freq,meanlowFreqOccips,errorlowFreqOccips)
    xlim([min(stimulation_frequencies)-2, max(stimulation_frequencies)+2]);
    title(['RLS Fit Results: ' parName ', run ' num2str(analRunNum) ' , low frequency report time points'])
    
    subplot(1,3,3);
    mseb(freqs_rls.freq,meanhighFreqOccips,errorhighFreqOccips)
    xlim([min(stimulation_frequencies)-2, max(stimulation_frequencies)+2]);
    title(['RLS Fit Results: ' parName ', run ' num2str(analRunNum) ' , high frequency report time points'])
    
    figure
    mseb(freqs_rls.freq,[meanlowFreqOccips; meanhighFreqOccips],[errorlowFreqOccips; errorhighFreqOccips])
    xlim([min(stimulation_frequencies)-2, max(stimulation_frequencies)+2]);
    
    
    
    end
end


  clearvars -except *frequencies  analRunNum  filenames  trialNum  trial_dur parName legend
end
