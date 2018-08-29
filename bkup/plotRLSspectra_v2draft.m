function plotRLSspectra_v2draft(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, trialNum, trial_dur,parName,legend,omniData)

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

%% Set up concatenated data

    % Find and remove artifacts
%     clean_data = removeArtifacts(cfg_trldef,'no');

    % Concatenate trials
%     concatData = concatenateTrials(all_data);
    
    % Create data files corresponding to perceptual reports
%     [lowFreqReportData, highFreqReportData] = epochOfInt('yes',omniData, concatData, pressLatency);
    
%% Set up nonConcat

    % Create data files corresponding to perceptual reports, keeping trials
    % seperate
    [lowFreqReportData, highFreqReportData] = epochOfInt('no',omniData, all_data, pressLatency);


    for latency = pressLatency; % 0:0.2:1.6; loop through possible button lags to look for effect
    loopCount = 1;
    for trial = 1:length(all_data.trial)
        for focus_electrode = analElectrodes %occipitals
            for dataSet = 1:3 %loop through each of the concatenated data sets (complete, low frequency report only, high frequency report only), results in a 3-plot subplot

            if dataSet == 1
                all_data = all_data; %concatData;
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

            trial_dur = ceil(max(all_data.time{1,trial}));
            
            %create structure for this trial
            trial_data = [];
            trial_data.label = all_data.label;
            trial_data.time{1,1} = all_data.time{1,trial};
            trial_data.trial{1,1} = all_data.trial{1,trial};
            trial_data.fsample = all_data.fsample;

            %% RLS Analysis
            rls_data = trial_data;
            single_rls(1) = trial_data;
            single_rls(2) = trial_data;
            msglength = 0;
            for iTrial = 1:size(trial_data.trial, 2);

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

            if dataSet == 1
                figure;
            end

            electrode_to_plot = focus_electrode; %focus_electrode
            for iTrial = 1
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
        
        allTimeTrial(trial,:) = mean(concatOccips);
        lowTimeTrial(trial,:) = mean(lowFreqOccips);
        highTimeTrial(trial,:) = mean(highFreqOccips);
        
        loopCount = 1;
        
        clear concatOccips lowFreqOccips highFreqOccips
        
    end
    
    % Plot the means with error
    
    errorallTimeTrial = ste(allTimeTrial);
    meanallTimeTrial = mean(allTimeTrial);
    
    errorlowTimeTrial = ste(lowTimeTrial);
    meanlowTimeTrial = mean(lowTimeTrial);
    
    errorhighTimeTrial = ste(highTimeTrial);
    meanhighTimeTrial = mean(highTimeTrial);
    
    figure
    subplot(1,3,1);
    mseb(freqs_rls.freq,meanallTimeTrial,errorallTimeTrial)
    xlim([min(stimulation_frequencies)-2, max(stimulation_frequencies)+2]);
    title(['RLS Fit Results: ' parName ', run ' num2str(analRunNum) ' , all time points'])
    
    subplot(1,3,2);
    mseb(freqs_rls.freq,meanlowTimeTrial,errorlowTimeTrial)
    xlim([min(stimulation_frequencies)-2, max(stimulation_frequencies)+2]);
    title(['RLS Fit Results: ' parName ', run ' num2str(analRunNum) ' , low frequency report time points'])
    
    subplot(1,3,3);
    mseb(freqs_rls.freq,meanhighTimeTrial,errorhighTimeTrial)
    xlim([min(stimulation_frequencies)-2, max(stimulation_frequencies)+2]);
    title(['RLS Fit Results: ' parName ', run ' num2str(analRunNum) ' , high frequency report time points'])
    
    figure
%     title(['Low vs. High Report: RLS, 
    mseb(freqs_rls.freq,[meanlowTimeTrial; meanhighTimeTrial],[errorlowTimeTrial; errorhighTimeTrial])
    xlim([min(stimulation_frequencies)-2, max(stimulation_frequencies)+2]);
%     legend('14 Hz Report','17 Hz Report');
    
    
    
    end


  clearvars -except *frequencies  analRunNum  filenames  trialNum  trial_dur parName legend
end
