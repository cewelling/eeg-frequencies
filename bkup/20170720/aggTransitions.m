function [ transitions ] = aggTransitions(parName, iPar, runName, date, EEGfile, paramsFlag)
% aggTransitions(parName, runName, date, EEGfile)  Saves chunks 
% of RLS data corresponding to participant-reported transitions.
%
% [transitions] = aggTransitions(parName, runName, date, EEGfile) returns a
% struct with entries corresponding to 'low to high', 'high to low', 'low
% to low', and 'high to high' transitions, respectively. 
%
% transitions struct fields: 
% type: e.g. high to low
% f1: low frequency traces
% f2: high frequency traces
% crit: matrix of criteria that can be used for exclusion
%
% Within the crit matrix, columns correspond to:
% 1: length of first dom state
% 2: length of second dom state
% 3: length of mixed state
% 4: length of largest gap
% 5: length of time between end of first and beginning of second dom state
% 6: SNR of the trial the transition occurred in
% 7: SNR of the low frequency
% 8: SNR of the high frequncy
% 9: timePoint of the transition
% 10: trial in which the transition occurred
%
% Plots average of transitions for the run designated. Note that
% transitions are not excluded on the basis of SNR at this point.
% 
% Called from: analysisController.m, makeParPlot.m
% Dependencies: analysisParams.m, runRLS.m, parseSwitchData.m

%% Set-up

% load parameters
analysisParams

% Transition codes to look for in list of epochs, and their labels
caseType = {'low to high'; 'high to low'; 'low to low'; 'high to high'};

% Preallocate space
caseTimes = nan(60,4);
caseTrials = nan(60,4);

%% Prepare to label transitions with trial-by-trial SNR values (for exclusion later)
 
% space to store SNRs for each trial
snrVect = nan(1, numTrials);
snrVect_lo = nan(1, numTrials);
snrVect_hi = nan(1, numTrials);

% ensure that SNRs have been calculated for this run
thisRunSNRs = dir(['pre-processing/highSNRelecs/' parName '_' runName '*' electrodeSet.name '.mat']);
if isempty(thisRunSNRs)
    runFFT(EEGfile, runName, parName, iPar, date);
end

% load and store SNRs
for iTrial = 1:numTrials
    if exist(['pre-processing/highSNRelecs/' parName '_' runName '_trial' num2str(iTrial) '_' electrodeSet.name '.mat'],'file')
        load(['pre-processing/highSNRelecs/' parName '_' runName '_trial' num2str(iTrial) '_' electrodeSet.name '.mat']);
        snrVect(iTrial) = maxSNRs(2, ismember(maxSNRs(1,:), snrElecs));
        snrVect_lo(iTrial) = lFreqSNRs(2, ismember(lFreqSNRs(1,:), snrElecs));
        snrVect_hi(iTrial) = hFreqSNRs(2, ismember(hFreqSNRs(1,:), snrElecs));
    else % trial doesn't exist (due to mistake)
        snrVect(iTrial) = NaN;
        snrVect_lo(iTrial) = NaN;
        snrVect_hi(iTrial) = NaN;
    end
end

%% Find and categorize transitions based on perceptual report

% use sim schedule for sim trials
if ~isempty(strfind(runName, 'Sim')) && useSimSchedule
    if ~exist(['transitions/tLists_simSched/' parName '_' runName '.mat'], 'file')
        if ~exist(['simSchedules/' parName '_' runName '.mat'], 'file')
            [~, simSchedule] = getEpochs(parName, runName, date, paramsFlag);
        else
            load(['simSchedules/' parName '_' runName '.mat'])
        end
        % find transitions using sim schedule
        [l2h, h2l, l2l, h2h] = findTransitions(simSchedule, parName, runName);
    else
        load(['transitions/tLists_simSched/' parName '_' runName '.mat'])
    end
    
else % use perceptual report
    if ~isempty(strfind(runName, 'Sim')) % get sim schedule for simulation runs
        if ~exist(['transitions/tLists_simSched/' parName '_' runName '.mat'], 'file')
            if ~exist(['simSchedules/' parName '_' runName '.mat'], 'file')
                [~, simSchedule] = getEpochs(parName, runName, date, paramsFlag);
            else
                load(['simSchedules/' parName '_' runName '.mat'])
            end
            % find transitions using sim schedule
            [l2h, h2l, l2l, h2h] = findTransitions(simSchedule, parName, runName);
        else
            load(['transitions/tLists_simSched/' parName '_' runName '.mat'])
        end
        
        % store sim schedule switches
        l2hSimSwitches = l2h{1};
        h2lSimSwitches = h2l{1};
        
    end
    if ~exist(['transitions/tLists_CERW/' parName '_' runName '.mat'], 'file')
        % find transitions using Caroline's switch data file
        [l2h, h2l, l2l, h2h] = parseSwitchData(parName, runName, paramsFlag);
    else
        load(['transitions/tLists_CERW/' parName '_' runName '.mat'])
    end
    
    % get simulation switches relative to button presses
    if ~isempty(strfind(runName, 'Sim'))
        l2hPressDelays = l2h{1} - l2hSimSwitches;
        h2lPressDelays = h2l{1} - h2lSimSwitches;
    end
    
end

%%% In l2h, h2l, etc.
% 1: transition points
% 2: transition trials
% 3: length of time between end of first and beginning of second dom state
% 4: length of first dom state
% 5: length of second dom state
% 6: length of mixed state 
% 7: length of largest gap

%%% Criteria %%%
% 1: length of first dom state
% 2: length of second dom state
% 3: length of mixed state
% 4: length of largest gap
% 5: length of time between end of first and beginning of second dom state
% 6: SNR of the trial the transition occurred in
% 7: SNR of the low frequency
% 8: SNR of the high frequency
% 9: transition time
% 10: transition trial

caseTimes(1:length(l2h{1}),1) = l2h{1};
caseTrials(1:length(l2h{2}),1) = l2h{2};
trialSNRs = snrVect(l2h{2}(:));
loSNRs = snrVect_lo(l2h{2}(:));
hiSNRs = snrVect_hi(l2h{2}(:));
caseCrit{1} = [l2h{4} l2h{5} l2h{6} l2h{7} l2h{3} trialSNRs' loSNRs' hiSNRs' l2h{1} l2h{2}];

caseTimes(1:length(h2l{1}),2) = h2l{1};
caseTrials(1:length(h2l{2}),2) = h2l{2};
trialSNRs = snrVect(h2l{2}(:));
loSNRs = snrVect_lo(h2l{2}(:));
hiSNRs = snrVect_hi(h2l{2}(:));
caseCrit{2} = [h2l{4} h2l{5} h2l{6} h2l{7} h2l{3} trialSNRs' loSNRs' hiSNRs' h2l{1} h2l{2}];

caseTimes(1:length(l2l{1}),3) = l2l{1};
caseTrials(1:length(l2l{2}),3) = l2l{2};
trialSNRs = snrVect(l2l{2}(:));
loSNRs = snrVect_lo(l2l{2}(:));
hiSNRs = snrVect_hi(l2l{2}(:));
caseCrit{3} = [l2l{4} l2l{5} l2l{6} l2l{7} l2l{3} trialSNRs' loSNRs' hiSNRs' l2l{1} l2l{2}];

caseTimes(1:length(h2h{1}),4) = h2h{1};
caseTrials(1:length(h2h{2}),4) = h2h{2};
trialSNRs = snrVect(h2h{2}(:));
loSNRs = snrVect_lo(h2h{2}(:));
hiSNRs = snrVect_hi(h2h{2}(:));
caseCrit{4} = [h2h{4} h2h{5} h2h{6} h2h{7} h2h{3} trialSNRs' loSNRs' hiSNRs' h2h{1} h2h{2}];

%% Get the smoothed RLS time-amplitude data for the times corresponding to the transitions (for averaged electrodes)


if exist(['rls_data/smoothedRLS/' parName '_' runName '_' analFreqLabel '_' electrodeSet.name '.mat'], 'file') 
    load(['rls_data/smoothedRLS/' parName '_' runName '_' analFreqLabel '_' electrodeSet.name '.mat'])
else
    [rls_data, rls_time] = runRLS(parName, runName, date, EEGfile);
end

%% Normalize the RLS data (for averaged electrodes)

normRLS_data = [];

% If selected, get baseline using average across adjacent frequencies 
if strcmp(normType, 'freqWin') && exist([blineDir parName '_' runName '_' blineLoc '.mat'], 'file')
    load([blineDir parName '_' runName '_' blineLoc '.mat']);
elseif strcmp(normType, 'freqWin')
    baselines = getBaseline(parName, runName);
end

for iTrial = 1:size(rls_data(1).amp, 2)
    
    for iFreq = 1:length(rls_data);
        
        RLSamp = (rls_data(iFreq).amp{iTrial})';
        
        % account for trials that weren't recorded properly (mistakes)
        if isempty(RLSamp)
            normRLSamp = [];
            continue;
        end
        
        % take normalizing/standardizing step
        if strcmp(normType, 'mean')
            normRLSamp = RLSamp - nanmean(RLSamp);
        elseif strcmp(normType, 'norm')
            normRLSamp = 2*(RLSamp - min(RLSamp))/(max(RLSamp) - min(RLSamp)) - 1;
        elseif strcmp(normType, 'z')
            normRLSamp = zscore(RLSamp);
        elseif strcmp(normType, 'none')
            normRLSamp = RLSamp;
        elseif strcmp(normType, 'noStimBase')
            normRLSamp = RLSamp - baselines(iFreq); %,iTrial);
        elseif strcmp(normType, 'freqWin')
            normRLSamp = RLSamp ./ nanmean(baselines{iTrial}{iFreq})';
        end
        
        normRLS_data(iFreq).amp{iTrial} = normRLSamp';
    end
end

%% Get pieces of RLS data that correspond to transitions

halfLength = sampRate*transHalf; 
fullLength = halfLength*2+1; % twice the half length plus the transition point

analElecs = electrodeSet.nums;

for iCase = 1:size(caseTrials,2) %for each type of transition
    
    thisCaseTrials = caseTrials(~isnan(caseTrials(:,iCase)),iCase); %list of trials containing this transition
    
    f1Epochs = cell(1,length(analElecs)); % low stim freq
    f2Epochs = cell(1,length(analElecs)); % high stim freq
    f1Epochs_elecAvg = [];
    f2Epochs_elecAvg = [];
    
    for iEpoch = 1:length(thisCaseTrials)
        
        transitionTime = caseTimes(iEpoch,iCase);
        
        for iElec = 0:length(analElecs)
            
            if iElec == 0 % average electrodes together
                
                RLSf1 = normRLS_data(1).amp{thisCaseTrials(iEpoch)};
                RLSf2 = normRLS_data(2).amp{thisCaseTrials(iEpoch)};
                
            else % handle individual electrodes
                
                %% Ensure that individual channel rls analysis has been completed for this run
                
                thisRunFiles = dir(['rls_data/perChannel/' parName '_' runName '*' num2str(analElecs(iElec)) '_' analFreqLabel '.mat']);
                if isempty(thisRunFiles)
                    error(['Run RLS analysis from analysisController with electrode ' num2str(analElecs(iElec)) ' and ' analFreqLabel ' selected']);
                end

                %% get smoothed RLS traces from appropriate trial
                if exist(['rls_data/perChannel/' parName '_' runName '_trial' num2str(thisCaseTrials(iEpoch)) '_freq1_' num2str(analElecs(iElec)) '_' analFreqLabel '.mat'], 'file')
                    load(['rls_data/perChannel/' parName '_' runName '_trial' num2str(thisCaseTrials(iEpoch)) '_freq1_' num2str(analElecs(iElec)) '_' analFreqLabel '.mat']);
                    RLSf1 = thisElecRLS;
                    load(['rls_data/perChannel/' parName '_' runName '_trial' num2str(thisCaseTrials(iEpoch)) '_freq2_' num2str(analElecs(iElec)) '_' analFreqLabel '.mat']);
                    RLSf2 = thisElecRLS;
                else
                    continue; % account for trials that weren't recorded properly (mistakes)
                end
                
                %% normalize the RLS traces
                if strcmp(normType, 'mean')
                    RLSf1 = RLSf1 - nanmean(RLSf1);
                    RLSf2 = RLSf2 - nanmean(RLSf2);
                elseif strcmp(normType, 'none')
                    % do nothing
                else
                    error(['Update code to handle ' normType ' normalization']);
                end
                
            end
            
            %% extract transition data
            
            % find the index of the transition time
            [~, tTimeIndex] = min(abs(rls_time - transitionTime));
            
            % epoch is cut off by beginning, do not include it (to avoid discontinuities in transition traces)
            if (tTimeIndex - halfLength) < 1
                epochRLSf1 = nan(1, 2*halfLength + 1);
                epochRLSf2 = nan(1, 2*halfLength + 1);
                % pad the beginning with NaNs
                %             epochRLSf1 = RLSf1(1:tTimeIndex + halfLength);
                %             epochRLSf2 = RLSf2(1:tTimeIndex + halfLength);
                %             padSize = fullLength - length(epochRLSf1);
                %
                %             epochRLSf1 = [nan(1,padSize), epochRLSf1];
                %             epochRLSf2 = [nan(1,padSize), epochRLSf2];
                
                % epoch is cut off by end, do not include it (to avoid discontinuities in transition traces)
            elseif (tTimeIndex + halfLength) > length(RLSf1) - (cutoffWindow / 1000 * 512); % remember end is padded due to shift forward
                epochRLSf1 = nan(1, 2*halfLength + 1);
                epochRLSf2 = nan(1, 2*halfLength + 1);
                % pad the end with NaNs
                %             epochRLSf1 = RLSf1(tTimeIndex - halfLength : length(RLSf1) - (cutoffWindow / 1000 * 512));
                %             epochRLSf2 = RLSf2(tTimeIndex - halfLength : length(RLSf2) - (cutoffWindow / 1000 * 512));
                %             padSize = fullLength - length(epochRLSf1);
                %
                %             epochRLSf1 = [epochRLSf1,nan(1,padSize)];
                %             epochRLSf2 = [epochRLSf2,nan(1,padSize)];
                
                % epoch is not cut off
            else
                epochRLSf1 = RLSf1(tTimeIndex-halfLength:tTimeIndex + halfLength);
                epochRLSf2 = RLSf2(tTimeIndex-halfLength:tTimeIndex + halfLength);
            end
            
            % Collect epochs
            if iElec == 0 % average electrodes
                f1Epochs_elecAvg = [f1Epochs_elecAvg; epochRLSf1];
                f2Epochs_elecAvg = [f2Epochs_elecAvg; epochRLSf2];
            else 
                f1Epochs{iElec} = [f1Epochs{iElec}; epochRLSf1];
                f2Epochs{iElec} = [f2Epochs{iElec}; epochRLSf2];
            end
            
        end
    end
    
    % store transitions with properties (i.e. dom length, gap length, etc.)
    for iElec = 0:length(analElecs)
        if iElec == 0 % average electrodes
            transitions_elecAvg(iCase).type = caseType(iCase);
            transitions_elecAvg(iCase).f1 = f1Epochs_elecAvg;
            transitions_elecAvg(iCase).f2 = f2Epochs_elecAvg;
            transitions_elecAvg(iCase).crit = caseCrit{iCase};
        else
            allElecTransitions{iElec}(iCase).type = caseType(iCase);
            allElecTransitions{iElec}(iCase).f1 = f1Epochs{iElec};
            allElecTransitions{iElec}(iCase).f2 = f2Epochs{iElec};
            allElecTransitions{iElec}(iCase).crit = caseCrit{iCase};
        end
    end
end

% save transitions
for iElec = 0:length(analElecs)
    if iElec == 0 % average electrodes
        transitions = transitions_elecAvg;
        if ~isempty(strfind(runName, 'Sim')) && useSimSchedule
            save(['transitions/tListsRLS_simSched/' parName '_' runName '_' normType '_' analFreqLabel '_' electrodeSet.name], 'transitions');
        else
            if ~isempty(strfind(runName, 'Sim'))
                save(['transitions/tListsRLS/' parName '_' runName '_' normType '_' analFreqLabel '_' electrodeSet.name], 'transitions', 'l2hPressDelays', 'h2lPressDelays');
            else
                save(['transitions/tListsRLS/' parName '_' runName '_' normType '_' analFreqLabel '_' electrodeSet.name], 'transitions');
            end
        end
        fprintf('Saving RLS traces of transitions with criteria\n')
    else
        thisElec = analElecs(iElec);
        transitions = allElecTransitions{iElec}; % transitions for this electrode
        if ~isempty(strfind(runName, 'Sim')) && useSimSchedule
            save(['transitions/tListsRLS_simSched/' parName '_' runName '_' normType '_' analFreqLabel '_' num2str(thisElec)], 'transitions');
        else
            save(['transitions/tListsRLS/' parName '_' runName '_' normType '_' analFreqLabel '_' num2str(thisElec)], 'transitions');
        end
        fprintf('Saving individual electrode RLS traces of transitions with criteria\n')
    end
end

%% Plotting

% plot and return averaged transitions
transitions = transitions_elecAvg;

if strcmp(transPlotOrNot, 'yes')
    for iCase = 1:size(transitions,2) % for each type of transition
        
        if isempty(transitions(iCase).f1)
            disp(['no instances of' caseType(iCase)])
            
        elseif size(caseTrials(~isnan(caseTrials(:,iCase)),iCase),1) < minInstances
            disp(['not enough instances of' caseType(iCase)])
            
        else
            meanF1trace = nanmean(transitions(iCase).f1,1);
            errorF1trace = ste(transitions(iCase).f1);
            
            meanF2trace = nanmean(transitions(iCase).f2,1);
            errorF2trace = ste(transitions(iCase).f2);
            
            figure
            title([parName ' ' runName ' ' caseType(iCase)])
            hold on
            mseb([(-transHalf):(1/512):transHalf],[meanF1trace; meanF2trace],[errorF1trace; errorF2trace],[],1);
            vline(0,'k','transition')
            xlabel('Time from button press (s)');
            ylabel('Amplitude');
            legend('Low Frequency', 'High Frequency');
            
        end
        
    end
end

end

