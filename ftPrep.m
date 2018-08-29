function [ all_data ] = ftPrep( parName, runName, EEGfile, numTrials, date, lowBound, highBound, paramsFlag,cfg_preproc)
%ftPrep prepares raw EEG data for FFT and RLS analysis.
% Returns pre-processed data with noisy electrodes eliminated (all_data)
% lowBound and highBound are the bounds of the band pass filter

% load parameters
analysisParams

%% Trial definition

[cfg_trldef] = ftDefineTrial(EEGfile, numTrials, trialDur);

%%  Identify noisy electrodes
if ~exist(['pre-processing/badElecs/' parName '_' runName '.mat'], 'file')
    
    % Pre-process data without re-referencing, with wide bandpass filter
    checkNoise_data = ftPreProc(cfg_trldef, parName, EEGfile, date, 0, 2, 70,cfg_preproc);
    % ftPreProc( parName, EEGfile, date, reref, lowBound, highBound)
    % Wide bandpass filter -- we want to eliminate bad electrodes for any
    % relevant frequency of interest
    
    % Remove trials that aren't real due to misplaced triggers
    checkNoise_data = ftRemoveFakeTriggers( checkNoise_data, EEGfile );
    
    % Noisy electrode detection
    [ badElecs ] = ftDetectBadElecs(checkNoise_data, runName);
    
    % Save bad electrodes
    save(['pre-processing/badElecs/' parName '_' runName], 'badElecs'); 
    
else
    load(['pre-processing/badElecs/' parName '_' runName '.mat']);
end

%% Pre-processing data again
% Note: currently, we remove bad electrodes AFTER re-referencing

% Pre-process data with re-referencing, with narrow bandpass filter
%           ftPreProc( parName, EEGfile, date, reref, lowBound, highBound)
all_data = ftPreProc(cfg_trldef, parName, EEGfile, date, 1, lowBound, highBound,cfg_preproc);

% Remove trials that aren't real due to misplaced triggers
all_data = ftRemoveFakeTriggers( all_data, EEGfile,numTrials );

%% Remove bad electrodes (on a trial-by-trial basis)
badCount = zeros(32,1);
for iTrial = 1:length(badElecs)
    for iElec = 1:length(badElecs{iTrial})
        badElec = badElecs{iTrial}(iElec).num;
        % replace bad electrode with nans
        all_data.trial{iTrial}(badElec,:) = nan(1,size(all_data.trial{iTrial}(badElec,:),2));
        
        badCount(badElec) = badCount(badElec) + 1;
    end
end

% if electrode was noisy in more than three trials, remove it altogether
for i = 1:length(badCount)
    if badCount(i) > 3
        for iTrial = 1:length(badElecs)
            all_data.trial{iTrial}(i,:) = nan(1,size(all_data.trial{iTrial}(i,:),2));
        end
    end
end

%% Test electrodes (for finding a metric for noisy electrodes)
elecsChosen = []; % include elecs to be analyzed, plotted
testElecs = {};
for iTrial = 1:length(all_data.trial)
    trialData = all_data.trial{iTrial};
    count = 1;
    for elec = elecsChosen;
        kur = kurtosis(trialData(elec,:));
        vr = var(trialData(elec,:));
        mn = mean(trialData(elec,:));
        testElecs{iTrial}(count).num = elec;
        testElecs{iTrial}(count).kur = kur;
        testElecs{iTrial}(count).var = vr;
        testElecs{iTrial}(count).mean = mn;
        count = count + 1;
    end
end
save(['pre-processing/testElecs/' parName '_' runName], 'testElecs');

end

