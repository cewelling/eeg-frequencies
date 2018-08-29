function [ rls_data, rls_time ] = runRLS(parName, runName, date, EEGfile)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

% Called by: rateFFT.m
% load parameters
analysisParams

%% Get stored analysis electrodes

rls_data = [];
msglength = 0;

analElecs = electrodeSet.nums;

for iFreq = 1:length(analFreqs) % size([stimFreqs imFreqs])
    
    %% Trial definition and pre-processing
    
    if exist(['pre-processing/preproc/' parName '_' runName '_freq' strrep(num2str(analFreqs(iFreq)),'.','-') '.mat'], 'file')
        load(['pre-processing/preproc/' parName '_' runName '_freq' strrep(num2str(analFreqs(iFreq)),'.','-') '.mat'])
    else
        try
            all_data = ftConfigParams(parName, runName, EEGfile, numTrials, trialDur, date, [], analFreqs(iFreq) - 0.5, analFreqs(iFreq) + 0.5);
        catch
            % Ensure that path is in correct order
            path(genpath('/Users/acspiegel/Documents/MATLAB/Toolboxes/FieldTrip'), path)
            all_data = ftConfigParams(parName, runName, EEGfile, numTrials, trialDur, date, [], analFreqs(iFreq) - 0.5, analFreqs(iFreq) + 0.5);
        end
        save(['pre-processing/preproc/' parName '_' runName '_freq' strrep(num2str(analFreqs(iFreq)),'.','-')], 'all_data');
    end
    
    for iTrial = 1:size(all_data.trial, 2); % # trials recorded in EEG (should = numTrials)
        
        fprintf(repmat('\b',1,msglength));
        msg = sprintf('Running RLS - Processing Trial No %d\n', iTrial);
        fprintf(msg);
        msglength = numel(msg);
        
        %for iFreq = 1:length(stimFreqs) % size([stimFreqs imFreqs])
        currFreq = analFreqs(iFreq);
        window = sampRate/currFreq*10; %1000/currFreq*5; (in timePoints)
        
        % Use the RLS filter
        [Hcn1,Hsn1] = AdaptiveRLSfast(window,(all_data.trial{iTrial}(analElecs,:))',currFreq,[],1,sampRate);
        
        cplx = complex(Hcn1,Hsn1);
        envelope = abs(cplx);
        envelope = envelope';
        
        for iElec = 1:length(analElecs)
            % smooth RLS amplitudes
            envelope(iElec, :) = HRsmoothing(envelope(iElec, :), 'gaussian', smooth_ksize, smooth_sd, 0);
            %envelope(iElec, :) = smoothts(envelope(iElec,:), 'g', smoothingWin,smoothingStd);
            %envelope(iElec, :) = smooth(envelope(iElec,:), smoothingPts);
            
            % eliminate first second of data, shift all remaining data over
            pad = round((cutoffWindow/1000)*sampRate);
            shiftedEnv(iElec, :) = [envelope(iElec, (pad+1):end) nan(1,pad)];
            %shiftedEnv(iElec, :) = [envelope(iElec, (pad+1):end) (nanmean(envelope((pad+1):end))*ones(1,pad))];
            % pads end with average of all the data
        end
        
        % Adjust trial numbering to account for mistakes
        thisTrial = iTrial;
        if ~isempty(strfind(EEGfile, 'cumulus05_marzRival1'))
            % account for loss of first two trials of this run
            thisTrial = thisTrial + 2;
        elseif ~isempty(strfind(EEGfile, 'stratus135_dartSim1'))
            % account for loss of first trial of this run
            thisTrial = thisTrial + 1;
        elseif ~isempty(strfind(EEGfile, 'cumulus01_dartSim3'))
            % account for loss of first trial of this run
            thisTrial = thisTrial +1;
        elseif ~isempty(strfind(EEGfile, 'cumulus01_marzRival3'))
            % account for loss of first two trials of this run
            thisTrial = thisTrial +2;
        end
        
        % save RLS traces for individual electrodes
        for iElec = 1:length(analElecs)
            thisElecRLS = shiftedEnv(iElec,:);
            save(['rls_data/perChannel/' parName '_' runName '_trial' num2str(thisTrial) ...
                '_freq' num2str(iFreq) '_' num2str(analElecs(iElec)) ...
                '_' analFreqLabel], 'thisElecRLS');
        end
        
        % average across electrodes
        shiftedEnv = nanmean(shiftedEnv, 1);
        rls_data(iFreq).amp {thisTrial} = shiftedEnv;
        rls_data(iFreq).electrodes {thisTrial} = electrodeSet.name;
    end
end

rls_time = all_data.time{1}; % time vector should be the same for all trials

% save RLS data (averaged across electrodes)
save(['rls_data/smoothedRLS/' parName '_' runName '_' analFreqLabel '_' electrodeSet.name], 'rls_data', 'rls_time');

%% Plot the time series, overlay perceptual epochs

if strcmp(RLSplotOrNot, 'yes')
    
    for iTrial = 1:size(all_data.trial, 2)
        
        trialTime = rls_time;
        
        % Average across electrodes
        lowBand = nanmean(rls_data(1).amp{iTrial}, 1);
        lowBandError = ste(rls_data(1).amp{iTrial});
        
        highBand = nanmean(rls_data(2).amp{iTrial}, 1);
        highBandError = ste(rls_data(2).amp{iTrial});
        
        % get list of participant's perceptual epochs
        [epochs, ~] = getEpochs(parName, runName, date);
        
        %plotting
        trialEpochs = epochs(epochs(:,1) == iTrial,:);
        
        figure
        title(['RLS: ', parName ',' runName ', trial ' num2str(iTrial)])
        hold on
        %         mseb(trialTime,[lowBand; highBand],[lowBandError; highBandError],[],1)
        plot(trialTime,lowBand,'b','linewidth',2)
        plot(trialTime,highBand,'r','linewidth',2)
        
        % find highest value that will be plotted
        checkMax = [max(lowBand+lowBandError) max(highBand+highBandError)];
        shadeLim = ceil(max(checkMax)); % round up
        
        for iEpoch = 1:size(trialEpochs,1)
            
            eStart = trialEpochs(iEpoch,2); % in seconds
            eEnd = trialEpochs(iEpoch,3); % in seconds
            horiz = [eStart eEnd];
            vert = [shadeLim shadeLim];
            
            %              startI = round(eStart*sampRate); % in data points
            %              endI = round(eEnd*sampRate); % in data points
            %              ePoints = endI - startI; % length of epoch in data points
            %              segFrac = 0.3; % proportion of epoch used for amplitude clustering
            
            % account for latency
            %              if startI > latency
            %                  startI_lat = startI - latency;
            %              else
            %                  startI_lat = 1;
            %              end
            
            %              segStart = startI_lat + round(ePoints*(1-segFrac)/2);
            %              seg = segStart:(segStart + round(ePoints*segFrac));
            
            if trialEpochs(iEpoch,5) > 0.5
                area(horiz,vert,'facecolor','r','facealpha',0.15,'edgealpha',0.1);
                %lOppAmps = [lOppAmps lowBand(seg)];
                %hCorrAmps = [hCorrAmps highBand(seg)];
            elseif trialEpochs(iEpoch,5) < -0.5
                area(horiz,vert,'facecolor','b','facealpha',0.15,'edgealpha',0.1);
                %lCorrAmps = [lCorrAmps lowBand(seg)];
                %hOppAmps = [hOppAmps highBand(seg)];
            else
                area(horiz,vert,'facecolor','y','facealpha',0.15,'edgealpha',0.1);
                %lMixedAmps = [lMixedAmps lowBand(seg)];
                %hMixedAmps = [hMixedAmps highBand(seg)];
            end
            
        end
        
    end
    
    %      % low frequency plot
    %      figure
    %      scatterX = [ones(1, length(lCorrAmps)) (ones(1, length(lMixedAmps)) + 1) (ones(1, length(lOppAmps)) + 2)];
    %      scatter(scatterX, [lCorrAmps lMixedAmps lOppAmps]);
end

end

