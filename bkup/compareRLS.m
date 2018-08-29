function compareRLS(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, trialNum, trial_dur,parName,legend,maxIndices,omniData)

%% Clear everything and establish where data is
% clearvars;
% clc;


%% Analysis Variables
for electrodeSet = 1 %1:2
    
    %% Load default anal params
    defaultAnalParams
    
    if electrodeSet == 1
        analElectrodes = electrodes;
        electrodeNames = 'allElectrodes';
    elseif electrodeSet == 2
        analElectrodes = occipitals;
        electrodeNames = 'occipitals';
    elseif electrodeSet == 3
        analElectrodes = maxIndices;
        electrodeNames = 'topResponders';
    end


    %% Practice Analysis with one subject
    fileID = fullfile(file_directory, filenames{1}(1).name);

    %% load config params
    [cfg_trldef, cfg_preproc] = configParams(fileID,discard_start,trial_dur,trialNum,analElectrodes)

    %% preprocess data
    all_data = ft_preprocessing(cfg_preproc);

    %% Percept Definition
    % Find all Percept Starts
    cfg_percdef.dataset = fileID;
    cfg_percdef.trialdef.eventtype = 'STATUS';
    cfg_percdef.trialfun = 'ft_trialfun_general';
    cfg_percdef.trialdef.prestim = 0;
    cfg_percdef.trialdef.poststim = 0;
    cfg_percdef.trialdef.eventvalue = 1:10;
    cfg_percdef = ft_definetrial(cfg_percdef);


%% Set up

    % Find and remove artifacts
%     clean_data = removeArtifacts(cfg_trldef,'yes');

%     % Concatenate trials
%     concatData = concatenateTrials(all_data);
%     
%     % Create data files corresponding to perceptual reports
%     [lowFreqReportData, highFreqReportData] = epochOfInt(omniData, concatData);
%     
%     all_data = lowFreqReportData;
%     trial_dur = ceil(max(all_data.time{1,1}));
%     trialNum = 1;


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

%             [Hcn1,Hsn1] = AdaptiveRLSfast(500,(all_data.trial{iTrial}(focus_electrode,:))',iFreq,[],1,512);
%             cplx = complex(Hcn1,Hsn1);
%             envelope = abs(cplx);
% 
%             rls_data.trial{iTrial} = envelope';
%             single_rls(j).amp{iTrial} = envelope';
        end


    end

    %% Plot the RLS Signal Amplitude Estimate
    
    trialDifferences = []; %store the difference scores
    meanTimes = mean([omniData(:,2) omniData(:,3)],2);
    meanTimes = [omniData(:,1), meanTimes];
    
    runLowDiff =  [];
    runHighDiff = [];
    
    for iTrial = 1:size(all_data.trial, 2);
%     for iTrial = [1 2 4 5 6]; %workaround for lack of button press data in 15-4

        for iFreq = 1:2;

            y{iFreq} = smooth(single_rls(iFreq).amp{iTrial}(focus_electrode, :), 150); 

        end
        
            thisOmni = omniData(omniData(:,1) == iTrial,:);

            if ~isempty(thisOmni)
                for iFreq = 1:2;

                    lowReportAmps = [];
                    highReportAmps = [];
                    for iTimePt = 1:length(single_rls(iFreq).time{iTrial})
                        [distance closeEpochRow] = min(abs(meanTimes(meanTimes(:,1) == iTrial,2)-single_rls(iFreq).time{iTrial}(1,iTimePt)));
                        if single_rls(iFreq).time{iTrial}(1,iTimePt) > thisOmni(closeEpochRow,2) && single_rls(iFreq).time{iTrial}(1,iTimePt) < thisOmni(closeEpochRow,3) && thisOmni(closeEpochRow,7) == 1
                            highReportAmps = [highReportAmps y{iFreq}(iTimePt)];
                        elseif single_rls(iFreq).time{iTrial}(1,iTimePt) > thisOmni(closeEpochRow,2) && single_rls(iFreq).time{iTrial}(1,iTimePt) < thisOmni(closeEpochRow,3) && thisOmni(closeEpochRow,7) == -1
                            lowReportAmps = [lowReportAmps y{iFreq}(iTimePt)];
                        end
                    end

                    if iFreq == 1 % if looking at the signal in the low frequency band, subtract the amplitude of the high from the low
                        runLowDiff = [runLowDiff (mean(lowReportAmps) - mean(highReportAmps))];
                    elseif iFreq == 2 % if looking at the signal in the low frequency band, subtract the amplitude of the low from the high
                        runHighDiff = [runHighDiff (mean(highReportAmps) - mean(lowReportAmps))];
                    end


                end
            end
        
        freqSeries{iTrial} = y;
            
        lowDiffs(iTrial,:) = y{1} - y{2};
        highDiffs(iTrial,:) = y{2} - y{1};
        

    end
    
    lowMinusHigh = mean(runLowDiff);
    lowMinusHighError = ste(runLowDiff);
    
    highMinusLow = mean(runHighDiff);
    highMinusLowError = ste(runHighDiff);
    
%     figure
%     title(['RLS Amplitude Differences: ', parName ', run ' num2str(analRunNum) ', focus electrode ' num2str(focus_electrode)])
%     hold on
%     bar(1:2,[lowMinusHigh highMinusLow])
%     Labels = {'14.166 Hz', '17 Hz'};
%     set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
%     errorbar(1:2,[lowMinusHigh highMinusLow],[lowMinusHighError highMinusLowError],'.')


    
    trialTime = single_rls(1).time{1}; % create variable for time, for plotting later

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

%     figure;
%     electrode_to_plot = focus_electrode; %focus_electrode
%     for iTrial = 1:trialNum
%         subplot(round(trialNum/4), 4, iTrial); hold on; % changed round to ceil: round(trialNum/4), 5
% 
%         plot(freqs_raw.freq, squeeze( freqs_raw.powspctrm(iTrial, electrode_to_plot, :) ), 'Color', colours{1},'LineWidth',2);
% 
%         plot(freqs_rls.freq, squeeze( freqs_rls.powspctrm(iTrial, electrode_to_plot, :) ), '--r');
%         xlim([min(stimulation_frequencies)-5, max(stimulation_frequencies)+5]);
%     %     ylim([0, 1]);
%     %     xlim([0, max([10, stimulation_frequencies])+10]);
%     end
% 
%     box off
%     title(['RLS Fit Results: ' parName ', run ' num2str(analRunNum) ' , ' electrodeNames ', ' num2str(trialNum) ' trials'])
%     figName = [plotDir parName '_run' num2str(analRunNum) '_RLSfits_' electrodeNames '.jpg']
%     saveas(gcf,figName,'jpg')

    
    %% SNR

    for iTrial = 1:size(all_data.trial, 2);
        for iFreq = 1:length(stimulation_frequencies)

%             snr{iTrial}(iFreq, 1:32) = ...
%                 freqs_raw.powspctrm(iTrial, :, find(freqs_raw.freq>stimulation_frequencies(iFreq), 1))...
%                 ./median(freqs_raw.powspctrm(iTrial, :, (freqs_raw.freq>stimulation_frequencies(iFreq)-3 & freqs_raw.freq<stimulation_frequencies(iFreq)+3)), 3);
%             snr{iTrial}(iFreq, 1:32) = ...
%                 freqs_raw.powspctrm(iTrial, :, find(freqs_raw.freq>stimulation_frequencies(iFreq), 1))...
%                 ./median(freqs_raw.powspctrm(iTrial, :, (freqs_raw.freq>stimulation_frequencies(iFreq)-0.167 & freqs_raw.freq<stimulation_frequencies(iFreq)+0.167)), 3);
            snr{iTrial}(iFreq, 1:32) = ...
                freqs_raw.powspctrm(iTrial, :, find(freqs_raw.freq>stimulation_frequencies(iFreq), 1))...
                ./median(freqs_raw.powspctrm(iTrial, :, (freqs_raw.freq>stimulation_frequencies(iFreq)-1 & freqs_raw.freq<stimulation_frequencies(iFreq)+1)), 3);

        end
    end
    
    SNRacrossTrials = [];
    
    for iTrial = 1:size(snr,2)
        trialSNR = mean(cell2mat(snr(iTrial)));
        SNRacrossTrials = [SNRacrossTrials; trialSNR];
    end
    
    SNRacrossTrials = mean(SNRacrossTrials);
        
    
%     [response topSNRelec] = max(mean(snr));
%     
%     topSNRelec
    
    % electrode group selection
    
    %option1
%     [SNRs,topResponders]=sort(SNRacrossTrials);
%     topResponders = topResponders((length(topResponders)-3):end)
    
    %option2
%     [SNRs,topResponders]=sort(snr,2);
%     topResponders = [topResponders(1,(size(topResponders,2)-1):size(topResponders,2)) topResponders(2,(size(topResponders,2)-1):size(topResponders,2))]
    
    %option3
    topResponders = [20 27 28 29 30];
%     topResponders = 23;

    toPlot = SNRacrossTrials;

    figure
    title('Max SNR Electrodes')
    hold on
    for iElec = 1:length(toPlot)
        h = bar(iElec,toPlot(iElec));
        if ismember(iElec,topResponders)
            set(h,'FaceColor','r')
        else
            set(h,'FaceColor','k')
        end
    end
    hold off


    % return;

end

% clearvars -except *frequencies  analRunNum  filenames  trialNum  trial_dur parName legend topResponders file_directory filenames

clearvars cfg* all_data freqs* legend rls_data single_rls y runHighDiff runLowDiff


%% Screen trials based on disproportional frequency response

for focus = topResponders(end)
    
    fileID = fullfile(file_directory, filenames{1}(1).name);

    %% load config params
    [cfg_trldef, cfg_preproc] = configParams(fileID,discard_start,trial_dur,trialNum,analElectrodes)

    %% preprocess data
    all_data = ft_preprocessing(cfg_preproc);

    %% Percept Definition
    % Find all Percept Starts
    cfg_percdef.dataset = fileID;
    cfg_percdef.trialdef.eventtype = 'STATUS';
    cfg_percdef.trialfun = 'ft_trialfun_general';
    cfg_percdef.trialdef.prestim = 0;
    cfg_percdef.trialdef.poststim = 0;
    cfg_percdef.trialdef.eventvalue = 1:10;
    cfg_percdef = ft_definetrial(cfg_percdef);


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

%         colours = {[83 148 255]/255, [255 117 117]/255};
        rls_data.trial{iTrial} = zeros(size(rls_data.trial{iTrial}));

        % Fundamental Stimulation Frequency
        for iFreq = stimulation_frequencies
            j = find(iFreq==stimulation_frequencies);
            cfg_rls.n_cycles = (trial_dur - discard_start) * iFreq;
            cfg_rls.stim_freq = iFreq;
            cfg_rls.channel = focus;
% 
            [single_rls(j).trial{iTrial}, single_rls(j).amp{iTrial}] = rls_slave( cfg_rls, all_data.trial{iTrial} );

            rls_data.trial{iTrial} = rls_data.trial{iTrial} + single_rls(j).trial{iTrial};
            
%             [Hcn1,Hsn1] = AdaptiveRLSfast(500,(all_data.trial{iTrial}(focus,:))',iFreq,[],1,512);
%             cplx = complex(Hcn1,Hsn1);
%             envelope = abs(cplx);
% 
%             rls_data.trial{iTrial} = envelope';
%             single_rls(j).amp{iTrial} = envelope';
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
    cfg_fft.channel = 1:32;

    freqs_raw = ft_freqanalysis(cfg_fft, all_data);
    freqs_rls = ft_freqanalysis(cfg_fft, rls_data);

    figure;
    electrode_to_plot = focus; %focus
    for iTrial = 1:trialNum
        subplot(round(trialNum/4), 4, iTrial); hold on; % changed round to ceil: round(trialNum/4), 5

        plot(freqs_raw.freq, squeeze( freqs_raw.powspctrm(iTrial, electrode_to_plot, :) ), 'Color', colours{1},'LineWidth',2);

        plot(freqs_rls.freq, squeeze( freqs_rls.powspctrm(iTrial, electrode_to_plot, :) ), '--r');
        xlim([min(stimulation_frequencies)-5, max(stimulation_frequencies)+5]);
    %     ylim([0, 1]);
    %     xlim([0, max([10, stimulation_frequencies])+10]);
    end
    suptitle(['FFT by trial: ' parName ', run ' num2str(analRunNum) ' , focus electrode: ' num2str(focus) ', ' num2str(trialNum) ' trials'])
    box off
%     figName = [plotDir parName '_run' num2str(analRunNum) '_RLSfits_' electrodeNames '.jpg']
%     saveas(gcf,figName,'jpg')
    
    % Determine ratio of frequency amplitude
    
    evenResponse = [];
    
    for iTrial = 1:trialNum
        
        responseRatio = freqs_raw.powspctrm(iTrial, focus, 426)/freqs_raw.powspctrm(iTrial, focus, 511);
        
        if responseRatio > 0.5 && responseRatio < 2
            evenResponse = [evenResponse iTrial];
        end
    end
        
        
    clearvars cfg* all_data freqs* legend rls_data single_rls y runHighDiff runLowDiff



end



%% Run through the same analysis for the top four SNR electrodes

multiElecLow = [];
multiElecHigh = [];

multiFreqSeries = [];

for focus = topResponders
    
    fileID = fullfile(file_directory, filenames{1}(1).name);

    %% load config params
    [cfg_trldef, cfg_preproc] = configParams(fileID,discard_start,trial_dur,trialNum,analElectrodes)

    %% preprocess data
    all_data = ft_preprocessing(cfg_preproc);

    %% Percept Definition
    % Find all Percept Starts
    cfg_percdef.dataset = fileID;
    cfg_percdef.trialdef.eventtype = 'STATUS';
    cfg_percdef.trialfun = 'ft_trialfun_general';
    cfg_percdef.trialdef.prestim = 0;
    cfg_percdef.trialdef.poststim = 0;
    cfg_percdef.trialdef.eventvalue = 1:10;
    cfg_percdef = ft_definetrial(cfg_percdef);


%% Set up

    % Find and remove artifacts
%     clean_data = removeArtifacts(cfg_trldef,'yes');

%     % Concatenate trials
%     concatData = concatenateTrials(all_data);
%     
%     % Create data files corresponding to perceptual reports
%     [lowFreqReportData, highFreqReportData] = epochOfInt(omniData, concatData);
%     
%     all_data = lowFreqReportData;
%     trial_dur = ceil(max(all_data.time{1,1}));
%     trialNum = 1;


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

%         colours = {[83 148 255]/255, [255 117 117]/255};
        rls_data.trial{iTrial} = zeros(size(rls_data.trial{iTrial}));

        % Fundamental Stimulation Frequency
        for iFreq = stimulation_frequencies
            j = find(iFreq==stimulation_frequencies);
            cfg_rls.n_cycles = (trial_dur - discard_start) * iFreq;
            cfg_rls.stim_freq = iFreq;
            cfg_rls.channel = focus;

            [single_rls(j).trial{iTrial}, single_rls(j).amp{iTrial}] = rls_slave( cfg_rls, all_data.trial{iTrial} );

            rls_data.trial{iTrial} = rls_data.trial{iTrial} + single_rls(j).trial{iTrial};
            
%             [Hcn1,Hsn1] = AdaptiveRLSfast(500,(all_data.trial{iTrial}(focus,:))',iFreq,[],1,512);
%             cplx = complex(Hcn1,Hsn1);
%             envelope = abs(cplx);
% 
%             rls_data.trial{iTrial} = envelope';
%             single_rls(j).amp{iTrial} = envelope';
        end


    end

    %% Plot the RLS Signal Amplitude Estimate

%     timeFig = figure;
    % scatterFig = figure;
    
    trialDifferences = []; %store the difference scores
    meanTimes = mean([omniData(:,2) omniData(:,3)],2);
    meanTimes = [omniData(:,1), meanTimes];
    
    runLowDiff =  [];
    runHighDiff = [];
    
    for iTrial = 1:size(all_data.trial, 2);
%     for iTrial = [1 2 4 5 6]; %workaround for lack of button press data in 15-4


%         figure(timeFig);
%         subplot(round(trialNum/2), 3, iTrial); hold on;

        for iFreq = 1:2;

            y{iFreq} = smooth(single_rls(iFreq).amp{iTrial}(focus, :), 150);        
%             plot(single_rls(iFreq).time{iTrial},y{iFreq}, 'Color', colours{iFreq});

        end
        
%         if ismember(iTrial,evenResponse) %accrue the means for only good trials
        if ismember(iTrial,[2 4]) %accrue the means for only good trials
        
            %for each of the time points in the rls estimate, determine if it
            %falls within the epoch report of high percpet
        
            thisOmni = omniData(omniData(:,1) == iTrial,:);

            if ~isempty(thisOmni)
                for iFreq = 1:2;

                    lowReportAmps = [];
                    highReportAmps = [];
                    for iTimePt = 1:length(single_rls(iFreq).time{iTrial})
                        [distance closeEpochRow] = min(abs(meanTimes(meanTimes(:,1) == iTrial,2)-single_rls(iFreq).time{iTrial}(1,iTimePt)));
                        if single_rls(iFreq).time{iTrial}(1,iTimePt) > thisOmni(closeEpochRow,2) && single_rls(iFreq).time{iTrial}(1,iTimePt) < thisOmni(closeEpochRow,3) && thisOmni(closeEpochRow,7) == 1
                            highReportAmps = [highReportAmps y{iFreq}(iTimePt)];
                        elseif single_rls(iFreq).time{iTrial}(1,iTimePt) > thisOmni(closeEpochRow,2) && single_rls(iFreq).time{iTrial}(1,iTimePt) < thisOmni(closeEpochRow,3) && thisOmni(closeEpochRow,7) == -1
                            lowReportAmps = [lowReportAmps y{iFreq}(iTimePt)];
                        end
                    end

                    if iFreq == 1 % if looking at the signal in the low frequency band, subtract the amplitude of the high from the low
                        runLowDiff = [runLowDiff (mean(lowReportAmps) - mean(highReportAmps))];
                    elseif iFreq == 2 % if looking at the signal in the low frequency band, subtract the amplitude of the low from the high
                        runHighDiff = [runHighDiff (mean(highReportAmps) - mean(lowReportAmps))];
                    end


                end
            end
        elseif ismember(iTrial,[1 3 5 6]) %accrue the means for only good trials
        
            %for each of the time points in the rls estimate, determine if it
            %falls within the epoch report of high percpet
        
            thisOmni = omniData(omniData(:,1) == iTrial,:);

            if ~isempty(thisOmni)
                for iFreq = 1:2;

                    lowReportAmps = [];
                    highReportAmps = [];
                    for iTimePt = 1:length(single_rls(iFreq).time{iTrial})
                        [distance closeEpochRow] = min(abs(meanTimes(meanTimes(:,1) == iTrial,2)-single_rls(iFreq).time{iTrial}(1,iTimePt)));
                        if single_rls(iFreq).time{iTrial}(1,iTimePt) > thisOmni(closeEpochRow,2) && single_rls(iFreq).time{iTrial}(1,iTimePt) < thisOmni(closeEpochRow,3) && thisOmni(closeEpochRow,7) == -1
                            highReportAmps = [highReportAmps y{iFreq}(iTimePt)];
                        elseif single_rls(iFreq).time{iTrial}(1,iTimePt) > thisOmni(closeEpochRow,2) && single_rls(iFreq).time{iTrial}(1,iTimePt) < thisOmni(closeEpochRow,3) && thisOmni(closeEpochRow,7) == 1
                            lowReportAmps = [lowReportAmps y{iFreq}(iTimePt)];
                        end
                    end

                    if iFreq == 1 % if looking at the signal in the low frequency band, subtract the amplitude of the high from the low
                        runLowDiff = [runLowDiff (mean(lowReportAmps) - mean(highReportAmps))];
                    elseif iFreq == 2 % if looking at the signal in the low frequency band, subtract the amplitude of the low from the high
                        runHighDiff = [runHighDiff (mean(highReportAmps) - mean(lowReportAmps))];
                    end


                end
            end
        end
        
        freqSeries{iTrial} = y;
            
        lowDiffs(iTrial,:) = y{1} - y{2};
        highDiffs(iTrial,:) = y{2} - y{1};
        

    end
    
    multiFreqSeries = [multiFreqSeries; freqSeries];
    
    lowMinusHigh = mean(runLowDiff);
    highMinusLow = mean(runHighDiff);
    
    multiElecLow = [multiElecLow; lowMinusHigh];
    multiElecHigh = [multiElecHigh; highMinusLow];
    
%     lowMinusHighError = ste(runLowDiff);
%     highMinusLowError = ste(runHighDiff);
%     
%     % make this a 3-panel bar plot, with three different portions of the
%     % button report section
%     
%     figure
%     title(['RLS Amplitude Differences: ', parName ', run ' num2str(analRunNum) ', focus electrode ' num2str(focus)])
%     hold on
%     bar(1:2,[lowMinusHigh highMinusLow])
%     Labels = {'14.166 Hz', '17 Hz'};
%     set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
%     errorbar(1:2,[lowMinusHigh highMinusLow],[lowMinusHighError highMinusLowError],'.')

    
    %% Plot difference / button correspondance
        
%     for iDiff = 1:size(all_data.trial, 2);
%         
%         trialEpochs = omniData(omniData(:,1) == iDiff,:);
%         
%         figure
%         title(['RLS time-frequency estimate: ', parName ', run ' num2str(analRunNum) ', trial ' num2str(iDiff)])
%         hold on
% %         legend('14 Hz Signal','17 Hz Signal','14 Hz Report','17 Hz Report')
%         plot(single_rls(1).time{1}, cell2mat(freqSeries{1,iDiff}(1)),'b','linewidth',2)
%         plot(single_rls(1).time{1}, cell2mat(freqSeries{1,iDiff}(2)),'r','linewidth',2)
%         
%         checkMax = [max(cell2mat(freqSeries{1,iDiff}(1))) max(cell2mat(freqSeries{1,iDiff}(2)))];
%         shadeLim = ceil(max(checkMax));
%         
%         for iEpoch = 1:size(trialEpochs,1)
%             
%             horiz = trialEpochs(iEpoch,2):trialEpochs(iEpoch,4):trialEpochs(iEpoch,3);
%             vert = [shadeLim shadeLim];
%             
%             if trialEpochs(iEpoch,7) > 0.5
%                 H = area(horiz,vert,'facecolor','r','facealpha',0.25,'edgealpha',0.1);
%             elseif trialEpochs(iEpoch,7) < -0.5
%                 H = area(horiz,vert,'facecolor','b','facealpha',0.25,'edgealpha',0.1);
%             else
%                 H = area(horiz,vert,'facecolor','y','facealpha',0.25,'edgealpha',0.1);
%             end
%         end
%     end
    
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
    electrode_to_plot = focus; %focus
    for iTrial = 1:trialNum
        subplot(round(trialNum/4), 4, iTrial); hold on; % changed round to ceil: round(trialNum/4), 5

        plot(freqs_raw.freq, squeeze( freqs_raw.powspctrm(iTrial, electrode_to_plot, :) ), 'Color', colours{1},'LineWidth',2);

        plot(freqs_rls.freq, squeeze( freqs_rls.powspctrm(iTrial, electrode_to_plot, :) ), '--r');
        xlim([min(stimulation_frequencies)-5, max(stimulation_frequencies)+5]);
    %     ylim([0, 1]);
    %     xlim([0, max([10, stimulation_frequencies])+10]);
    end

    box off
    title(['RLS Fit Results: ' parName ', run ' num2str(analRunNum) ' , ' electrodeNames ', ' num2str(trialNum) ' trials'])
    figName = [plotDir parName '_run' num2str(analRunNum) '_RLSfits_' electrodeNames '.jpg']
    saveas(gcf,figName,'jpg')
    
    
    clearvars cfg* all_data freqs* legend rls_data single_rls y runHighDiff runLowDiff



end
    
    %% Amplitude Difference

    lowMinusHigh = mean(multiElecLow);
    highMinusLow = mean(multiElecHigh);
    
    lowMinusHighError = ste(multiElecLow);
    highMinusLowError = ste(multiElecHigh);
    
    figure
%     title(['RLS Amplitude Differences: ', parName ', run ' num2str(analRunNum) ', focus electrodes: ' num2str(topResponders) ', trials included: ' num2str(evenResponse)])
    title(['RLS Amplitude Differences: ', parName ', run ' num2str(analRunNum) ', focus electrodes: ' num2str(topResponders)])
    hold on
    bar(1:2,[lowMinusHigh highMinusLow])
    Labels = {'14.166 Hz', '17 Hz'};
    set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
    errorbar(1:2,[lowMinusHigh highMinusLow],[lowMinusHighError highMinusLowError],'.')
    
    %% Time Series
    
    lowFreq = [];
    highFreq = [];
    
    lowLineProps.col = 'b';
    highLineProps.col = 'r';
    
    
    
    for iTrial = 1:size(multiFreqSeries,2)
%     for iTrial = [1 2 4 5 6]; %workaround for lack of button press data in 15-4

        
        for iElec = 1:size(multiFreqSeries,1)
            thisElecLow = cell2mat(multiFreqSeries{iElec,iTrial}(1));
            thisElecHigh = cell2mat(multiFreqSeries{iElec,iTrial}(2));
            
            lowFreq = [lowFreq; thisElecLow'];
            highFreq = [highFreq; thisElecHigh'];
        end
        
        lowBand = mean(lowFreq,1);
        lowBandError = ste(lowFreq);
        
        highBand = mean(highFreq,1);
        highBandError = ste(highFreq);
        
        %plotting
        
        trialEpochs = omniData(omniData(:,1) == iTrial,:);
        
        figure
        title(['RLS: ', parName ', run ' num2str(analRunNum) ', trial ' num2str(iTrial)])
        hold on
        mseb(trialTime,[lowBand; highBand],[lowBandError; highBandError],[],1)
%         plot(trialTime,lowBand,'b','linewidth',2)
%         plot(trialTime,highBand,'r','linewidth',2)
        
        checkMax = [max(lowBand+lowBandError) max(highBand+highBandError)];
        shadeLim = ceil(max(checkMax));
        
        for iEpoch = 1:size(trialEpochs,1)
            
            horiz = trialEpochs(iEpoch,2):trialEpochs(iEpoch,4):trialEpochs(iEpoch,3);
            vert = [shadeLim shadeLim];
            
            if trialEpochs(iEpoch,7) > 0.5
                H = area(horiz,vert,'facecolor','r','facealpha',0.15,'edgealpha',0.1);
            elseif trialEpochs(iEpoch,7) < -0.5
                H = area(horiz,vert,'facecolor','b','facealpha',0.15,'edgealpha',0.1);
            else
                H = area(horiz,vert,'facecolor','y','facealpha',0.15,'edgealpha',0.1);
            end
        end
        
    end
            
            
        

end
    %% TFR Analysis
    % Trial Definition
%     cfg_tfr = [];
%     cfg_preproc = [];
%     discard_start = 0;
%     [cfg_tfr, cfg_preproc] = configParams(fileID,discard_start,trial_dur,trialNum,analElectrodes);
%     all_data = ft_preprocessing(cfg_preproc);
% 
%     cfg_tfr = [];
%     cfg_tfr.channel = analElectrodes; %'all' or analElectrodes
%     cfg_tfr.foi = 0:(1/30):maxFreq; %frequency of interest  %original: step by 0.25
%     cfg_tfr.toi = -3:slidingWindowStep:trial_dur; %time of interest
%     cfg_tfr.keeptrials = 'yes'; % this makes sure all trials are separate
% 
%     cfg_tfr.method = 'mtmconvol'; % cfg_tfr.method = 'wavelet';
%     cfg_tfr.taper = 'hanning'; %sliding window shape (tapered)
%     cfg_tfr.t_ftimwin = slidingWindowWidth*ones(size(cfg_tfr.foi)); % cfg_tfr.width = 16;
% 
%     tfr_raw = ft_freqanalysis(cfg_tfr, all_data);
%     tfr_rls = ft_freqanalysis(cfg_tfr, rls_data);
%     
%     %% Save RLS data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     dataNameRAW = [rlsDir parName '_RLSRAW_run' num2str(analRunNum) '_' electrodeNames '.mat'];
%     dataNameRLS = [rlsDir parName '_RLSModel_run' num2str(analRunNum) '_' electrodeNames '.mat'];
%     
%     save(dataNameRAW,'tfr_raw')
%     save(dataNameRLS,'tfr_rls')
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%
%     
%     % for iFreq = 1:2
%     %  tfr_single_freq{iFreq} = ft_freqanalysis(cfg_tfr, single_rls(iFreq));
%     % end
%     cmap = cbrewer('seq', 'PuBu', 200);
% 
%     yLims = [1 0.8; 1 1; 1 1]; %trial 3
%     zLims = [0.01 .1; 0.01 .07; 0.01 0.07]; %trial 5
%     for iTrial = 1:trialNum
%         figure;
%         for iFreq = 1:length(stimulation_frequencies)
%             subplot(2, 1, iFreq);
%             cfg_plot = [];
%             cfg_plot.channel = freqs_raw.label(focus_electrode);
%             cfg_plot.channel = 'all';
%             cfg_plot.interactive = 'no';
%             cfg_plot.trials = iTrial;
%             cfg_plot.xlim = [0 trial_dur];
%             cfg_plot.ylim = [stimulation_frequencies(iFreq)-yLims(iFreq,1), stimulation_frequencies(iFreq)+yLims(iFreq,2)];
%             %cfg_plot.zlim = [zLims(iFreq,1),zLims(iFreq,2)];
%             %cfg_plot.baseline = [-3, -1];
%             cfg_plot.baselinetype = 'relative';
% 
%             gridxy2([], stimulation_frequencies);
% 
%             ft_singleplotTFR(cfg_plot, tfr_rls  );
%             %colormap(cmap);
%         end
%         box off
%         title(['RLS Heat Results: ' parName ', run ' num2str(analRunNum) ' , ' electrodeNames ', ' num2str(trialNum) ' trials'])
%         figName = [plotDir parName '_run' num2str(analRunNum) '_RLSheat_' electrodeNames '.jpg'];
%         saveas(gcf,figName,'jpg')
% 
%     end
% 
%   clearvars -except *frequencies  analRunNum  filenames  trialNum  trial_dur parName legend
% end
