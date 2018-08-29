function [multiFreqSeries, trialTime] = compareRLSfast(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, trialNum, trial_dur,parName,legend,omniData,timeWindow,analElect)

%% Clear everything and establish where data is
% clearvars;
% clc;

    defaultAnalParams
    analElectrodes = electrodes;
    
    focus_electrode = analElect; % set the analysis electrode

%     %% Practice Analysis with one subject
%     fileID = fullfile(file_directory, filenames{1}(1).name);
% 
%     %% load config params
%     [cfg_trldef, cfg_preproc] = configParams(fileID,discard_start,trial_dur,trialNum,analElectrodes)
% 
%     %% preprocess data
%     all_data = ft_preprocessing(cfg_preproc);
% 
%     %% Percept Definition
%     % Find all Percept Starts
%     cfg_percdef.dataset = fileID;
%     cfg_percdef.trialdef.eventtype = 'STATUS';
%     cfg_percdef.trialfun = 'ft_trialfun_general';
%     cfg_percdef.trialdef.prestim = 0;
%     cfg_percdef.trialdef.poststim = 0;
%     cfg_percdef.trialdef.eventvalue = 1:10;
%     cfg_percdef = ft_definetrial(cfg_percdef);
% 
% 
% %% Set up
% 
%     % Find and remove artifacts
% %     clean_data = removeArtifacts(cfg_trldef,'yes');
% 
% %     % Concatenate trials
% %     concatData = concatenateTrials(all_data);
% %     
% %     % Create data files corresponding to perceptual reports
% %     [lowFreqReportData, highFreqReportData] = epochOfInt(omniData, concatData);
% %     
% %     all_data = lowFreqReportData;
% %     trial_dur = ceil(max(all_data.time{1,1}));
% %     trialNum = 1;
% 
% 
%     %% RLS Analysis
%     rls_data = all_data;
%     single_rls(1) = all_data;
%     single_rls(2) = all_data;
%     msglength = 0;
%     for iTrial = 1:size(all_data.trial, 2);
% 
%         fprintf(repmat('\b',1,msglength));
%         msg = sprintf('Running RLS - Processing Trial No %d\n', iTrial);
%         fprintf(msg);
%         msglength = numel(msg);
% 
%         colours = {[83 148 255]/255, [255 117 117]/255};
%         rls_data.trial{iTrial} = zeros(size(rls_data.trial{iTrial}));
% 
%         % Fundamental Stimulation Frequency
%         for iFreq = stimulation_frequencies
%             j = find(iFreq==stimulation_frequencies);
%             cfg_rls.n_cycles = (trial_dur - discard_start) * iFreq;
%             cfg_rls.stim_freq = iFreq;
%             cfg_rls.channel = focus_electrode;
%             window = timeWindow; %msec %500
%             sampleRate = 512;
% 
%             [Hcn1,Hsn1] = AdaptiveRLSfast(window,(all_data.trial{iTrial}(focus_electrode,:))',iFreq,[],1,sampleRate);
%             cplx = complex(Hcn1,Hsn1);
%             envelope = abs(cplx);
%             envelope = envelope';
%             pad = (window/1000)*sampleRate;
%             
% %             timeShifted = [envelope((pad+1):end) zeros(1,pad)];
%             timeShifted = [envelope((pad+1):end) (mean(envelope((pad+1):end))*ones(1,pad))];
% 
%             single_rls(j).amp{iTrial} = timeShifted;
%         end
% 
% 
%     end
% 
%     %% Plot the RLS Signal Amplitude Estimate
%     
%     trialDifferences = []; %store the difference scores
%     meanTimes = mean([omniData(:,2) omniData(:,3)],2);
%     meanTimes = [omniData(:,1), meanTimes];
%     
%     runLowDiff =  [];
%     runHighDiff = [];
%     
%     for iTrial = 1:size(all_data.trial, 2);
% %     for iTrial = [1 2 4 5 6]; %workaround for lack of button press data in 15-4
% 
%         for iFreq = 1:2;
% 
%             y{iFreq} = smooth(single_rls(iFreq).amp{iTrial}, 150); 
% %             y{iFreq} = smoothts(single_rls(iFreq).amp{iTrial},'g',(2*window),window);
% 
%         end
%         
%             thisOmni = omniData(omniData(:,1) == iTrial,:);
% 
%             if ~isempty(thisOmni)
%                 for iFreq = 1:2;
% 
%                     lowReportAmps = [];
%                     highReportAmps = [];
%                     for iTimePt = 1:length(single_rls(iFreq).time{iTrial})
%                         [distance closeEpochRow] = min(abs(meanTimes(meanTimes(:,1) == iTrial,2)-single_rls(iFreq).time{iTrial}(1,iTimePt)));
%                         if single_rls(iFreq).time{iTrial}(1,iTimePt) > thisOmni(closeEpochRow,2) && single_rls(iFreq).time{iTrial}(1,iTimePt) < thisOmni(closeEpochRow,3) && thisOmni(closeEpochRow,7) == 1
%                             highReportAmps = [highReportAmps y{iFreq}(iTimePt)];
%                         elseif single_rls(iFreq).time{iTrial}(1,iTimePt) > thisOmni(closeEpochRow,2) && single_rls(iFreq).time{iTrial}(1,iTimePt) < thisOmni(closeEpochRow,3) && thisOmni(closeEpochRow,7) == -1
%                             lowReportAmps = [lowReportAmps y{iFreq}(iTimePt)];
%                         end
%                     end
% 
%                     if iFreq == 1 % if looking at the signal in the low frequency band, subtract the amplitude of the high from the low
%                         runLowDiff = [runLowDiff (mean(lowReportAmps) - mean(highReportAmps))];
%                     elseif iFreq == 2 % if looking at the signal in the low frequency band, subtract the amplitude of the low from the high
%                         runHighDiff = [runHighDiff (mean(highReportAmps) - mean(lowReportAmps))];
%                     end
% 
% 
%                 end
%             end
%         
%         freqSeries{iTrial} = y;
%             
%         lowDiffs(iTrial,:) = y{1} - y{2};
%         highDiffs(iTrial,:) = y{2} - y{1};
%         
% 
%     end
%     
%     lowMinusHigh = mean(runLowDiff);
%     lowMinusHighError = ste(runLowDiff);
%     
%     highMinusLow = mean(runHighDiff);
%     highMinusLowError = ste(runHighDiff);
%     
%     trialTime = single_rls(1).time{1}; % create variable for time, for plotting later
%     
% %     topResponders = [20 27 28 29 30];
%     topResponders = [29];
%     
%     clearvars cfg* all_data freqs* legend rls_data single_rls y runHighDiff runLowDiff
%     
%     
%     
%     
%     %% Run through the same analysis for the top four SNR electrodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

multiElecLow = [];
multiElecHigh = [];

multiFreqSeries = [];

for focus = focus_electrode %topResponders
    
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
            window = timeWindow; %msec %500 %I think this is 1000ms now -Alina
            sampleRate = 512;
            
            % Use the RLS filter
            [Hcn1,Hsn1] = AdaptiveRLSfast(window,(all_data.trial{iTrial}(focus_electrode,:))',iFreq,[],1,sampleRate);
            cplx = complex(Hcn1,Hsn1);
            envelope = abs(cplx);
            envelope = envelope';
            
            % Calculate SNR
%             % Compute noise
%             noiseFreqs = [iFreq - .0056 : .001 : iFreq + .0056]; % bandwidth identical to Zhang et. al., 2011
%             noiseAmps = zeros(length(noiseFreqs),length(envelope));
%             for i = [1 : 1 : length(noiseFreqs)]
%                 [Hcn1N, Hsn1N] = AdaptiveRLSfast(window,(all_data.trial{iTrial}(focus_electrode,:))',noiseFreqs(i),[],1,sampleRate);
%                 cplxN = complex(Hcn1N, Hsn1N);
%                 envelopeN = abs(cplxN);
%                 noiseAmps(i,:) = envelopeN;
%             end          
%             meanNoiseAmp = mean(noiseAmps, 1);
%             SNRs = envelope./meanNoiseAmp;
%             SNR = mean(SNRs);
%             SNRdev = std(SNRs);

%             SNRs = envelope./y_rms;
%             SNR = mean(SNRs);
%             SNRdev = std(SNRs);
            
            pad = (window/1000)*sampleRate;
            % Eliminates first second (window) of data, shifts data back, pads the end with the
            % average of all the data
%             timeShifted = [envelope((pad+1):end) zeros(1,pad)];
            timeShifted = [envelope((pad+1):end) (mean(envelope((pad+1):end))*ones(1,pad))];

            single_rls(j).amp{iTrial} = timeShifted;
            
        end


    end

    %% Plot the RLS Signal Amplitude Estimate
    
    trialDifferences = []; %store the difference scores
    meanTimes = mean([omniData(:,2) omniData(:,3)],2); % time in the middle of each perceptual epoch
    meanTimes = [omniData(:,1), meanTimes];
    % Reminder: omniData = [omniTrial, omniStart, omniStop, omniTime,
    % omniRow/3, omniKey, omniFreq];
    
    runLowDiff =  [];
    runHighDiff = [];
    
    for iTrial = 1:size(all_data.trial, 2);
%     for iTrial = [1 2 4 5 6]; %workaround for lack of button press data in 15-4


        for iFreq = 1:2;

            y{iFreq} = smooth(single_rls(iFreq).amp{iTrial}, 150); % uses 150 point moving average
%             y{iFreq} = smoothts(single_rls(iFreq).amp{iTrial},'g',(2*window),window);

        end
        
%         if ismember(iTrial,evenResponse) %accrue the means for only good trials
%         if ismember(iTrial,[2 4]) %accrue the means for only good trials
        
            %for each of the time points in the rls estimate, determine if it
            %falls within the epoch report of high percpet
        
            thisOmni = omniData(omniData(:,1) == iTrial,:);

            if ~isempty(thisOmni)
                for iFreq = 1:2;

                    lowReportAmps = [];
                    highReportAmps = [];
                    for iTimePt = 1:length(single_rls(iFreq).time{iTrial})
                        % Find the closest epoch to the current amplitude time point 
                        [distance closeEpochRow] = min(abs(meanTimes(meanTimes(:,1) == iTrial,2)-single_rls(iFreq).time{iTrial}(1,iTimePt)));
                        % Mark time point as either high or low frequency according to participant's perceptual report
                        if single_rls(iFreq).time{iTrial}(1,iTimePt) > thisOmni(closeEpochRow,2) && single_rls(iFreq).time{iTrial}(1,iTimePt) < thisOmni(closeEpochRow,3) && thisOmni(closeEpochRow,7) == 1
                            highReportAmps = [highReportAmps y{iFreq}(iTimePt)];
                        elseif single_rls(iFreq).time{iTrial}(1,iTimePt) > thisOmni(closeEpochRow,2) && single_rls(iFreq).time{iTrial}(1,iTimePt) < thisOmni(closeEpochRow,3) && thisOmni(closeEpochRow,7) == -1
                            lowReportAmps = [lowReportAmps y{iFreq}(iTimePt)];
                        end
                    end

                    if iFreq == 1 % if looking at the signal in the low frequency band, subtract the amplitude of the high from the low
                        runLowDiff = [runLowDiff (mean(lowReportAmps) - mean(highReportAmps))]; % We expect this to be positive 
                    elseif iFreq == 2 % if looking at the signal in the high frequency band, subtract the amplitude of the low from the high
                        runHighDiff = [runHighDiff (mean(highReportAmps) - mean(lowReportAmps))];
                    end


                end
            end
        
        freqSeries{iTrial} = y; % Reminder: y has the smoothed amplitude estimates
            
        lowDiffs(iTrial,:) = y{1} - y{2}; % difference in amplitudes of frequencies
        highDiffs(iTrial,:) = y{2} - y{1};
        

    end
    
    multiFreqSeries = [multiFreqSeries; freqSeries];
    
    lowMinusHigh = mean(runLowDiff);
    highMinusLow = mean(runHighDiff);
    
    multiElecLow = [multiElecLow; lowMinusHigh];
    multiElecHigh = [multiElecHigh; highMinusLow];
    
    trialTime = single_rls(1).time{1}; % create variable for time, for plotting later

    
    
    clearvars cfg* all_data freqs* legend rls_data single_rls y runHighDiff runLowDiff
end

    
    %% Amplitude Difference

    lowMinusHigh = mean(multiElecLow);
    highMinusLow = mean(multiElecHigh);
    
    lowMinusHighError = ste(multiElecLow);
    highMinusLowError = ste(multiElecHigh);
    
    figure
%     title(['RLS Amplitude Differences: ', parName ', run ' num2str(analRunNum) ', focus electrodes: ' num2str(topResponders) ', trials included: ' num2str(evenResponse)])
    title(['RLS Amplitude Differences: ', parName ', run ' num2str(analRunNum) ', focus electrodes: ' num2str(analElect)])
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
        
        % Would average across electrodes if there were multiple electrodes
        lowBand = mean(lowFreq,1);
        lowBandError = ste(lowFreq);
        
        highBand = mean(highFreq,1);
        highBandError = ste(highFreq);
        
        %plotting
        
        trialEpochs = omniData(omniData(:,1) == iTrial,:);
        
        figure
        title(['RLS: ', parName ', run ' num2str(analRunNum) ', trial ' num2str(iTrial)])
        hold on
%         mseb(trialTime,[lowBand; highBand],[lowBandError; highBandError],[],1)
        plot(trialTime,lowBand,'b','linewidth',2)
        plot(trialTime,highBand,'r','linewidth',2)
        
        % find highest value that will be plotted 
        checkMax = [max(lowBand+lowBandError) max(highBand+highBandError)];
        shadeLim = ceil(max(checkMax)); % round up
        
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