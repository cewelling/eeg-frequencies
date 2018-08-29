function [ transitions ] = aggTransitions(parName, runName, date, EEGfile)
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
caseCode = [-1 1; 1 -1; 0 1; 0 -1; 1 0; -1 0];
caseType = {'low to high'; 'high to low'; 'low to low'; 'high to high'};

% Preallocate space
caseTimes = zeros(50,4);
caseTrials = zeros(50,4);

%% Prepare to label transitions with trial-by-trial SNR values (for exclusion later)

% Get stored analysis electrodes 
 if strfind(runName, 'dart')
     load(['pre-processing/highSNRelecs/' parName '_darts_' num2str(commonElecs) 'of' num2str(numElecs) 'elecs']);
 else
     load(['pre-processing/highSNRelecs/' parName '_marz_' num2str(commonElecs) 'of' num2str(numElecs) 'elecs']);
 end
 
% load SNRs for each trial in this run
snrVect = zeros(1, numTrials);
for iTrial = 1:numTrials
    load(['pre-processing/highSNRelecs/' parName '_' runName '_trial' num2str(iTrial) '_' num2str(numElecs) 'elecs.mat']);
    snrVect(iTrial) = maxSNRs(2, ismember(maxSNRs(1,:),elecs));
end

%% Find and categorize transitions based on perceptual report

if ~exist(['transitions/tLists_CERW/' parName '_' runName '.mat'], 'file')
    
    % find transitions using Caroline's switch data file
    [l2h, h2l, l2l, h2h] = parseSwitchData(parName, runName);
else
    load(['transitions/tLists_CERW/' parName '_' runName '.mat'])
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

caseTimes(1:length(l2h{1}),1) = l2h{1};
caseTrials(1:length(l2h{2}),1) = l2h{2};
trialSNRs = snrVect(l2h{2}(:));
caseCrit{1} = [l2h{4} l2h{5} l2h{6} l2h{7} l2h{3} trialSNRs'];

caseTimes(1:length(h2l{1}),2) = h2l{1};
caseTrials(1:length(h2l{2}),2) = h2l{2};
trialSNRs = snrVect(h2l{2}(:));
caseCrit{2} = [h2l{4} h2l{5} h2l{6} h2l{7} h2l{3} trialSNRs'];

caseTimes(1:length(l2l{1}),3) = l2l{1};
caseTrials(1:length(l2l{2}),3) = l2l{2};
trialSNRs = snrVect(l2l{2}(:));
caseCrit{3} = [l2l{4} l2l{5} l2l{6} l2l{7} l2l{3} trialSNRs'];

caseTimes(1:length(h2h{1}),4) = h2h{1};
caseTrials(1:length(h2h{2}),4) = h2h{2};
trialSNRs = snrVect(h2h{2}(:));
caseCrit{4} = [h2h{4} h2h{5} h2h{6} h2h{7} h2h{3} trialSNRs'];

%% Get the smoothed RLS time-amplitude data for the times corresponding to the transitions


if exist(['rls_data/smoothedRLS/' parName '_' runName '.mat'], 'file') && isequal(analFreqs, stimFreqs)
    load(['rls_data/smoothedRLS/' parName '_' runName '.mat'])
elseif exist(['rls_data/smoothedRLS/' parName '_' runName '_harmonics.mat'], 'file') && isequal(analFreqs,harFreqs)
    load(['rls_data/smoothedRLS/' parName '_' runName '_harmonics.mat'])
else
    [rls_data, rls_time] = runRLS(parName, runName, date, EEGfile);
end

%% Normalize the RLS data

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
            normRLSamp = RLSamp - mean(RLSamp);
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

for iCase = 1:size(caseTrials,2) %for each type of transition
    
    thisCaseTrials = nonzeros(caseTrials(:,iCase)); %list of trials containing this transition
    
    f1Epochs = []; % low stim freq
    f2Epochs = []; % high stim freq
    
    for iEpoch = 1:length(thisCaseTrials)
        
        transitionTime = caseTimes(iEpoch,iCase);
        
        % get RLS data from appropriate trial
        RLSf1 = normRLS_data(1).amp{thisCaseTrials(iEpoch)};
        RLSf2 = normRLS_data(2).amp{thisCaseTrials(iEpoch)};
        
        %         if isempty(RLSf1)
        %             continue;
        %         end
        
        % find the index of the transition time
        [~, tTimeIndex] = min(abs(rls_time - transitionTime));
        
        % extract transition data
        
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
        f1Epochs = [f1Epochs; epochRLSf1];
        f2Epochs = [f2Epochs; epochRLSf2];
        
    end
    
    % store transitions with properties (i.e. dom length, gap length, etc.)
    transitions(iCase).type = caseType(iCase);
    transitions(iCase).f1 = f1Epochs;
    transitions(iCase).f2 = f2Epochs;
    transitions(iCase).crit = caseCrit{iCase};    
end

% save transitions
if isequal(analFreqs, stimFreqs)
    save(['transitions/tListsRLS/' parName '_' runName '_' normType], 'transitions');
elseif isequal(analFreqs, harFreqs)
    save(['transitions/tListsRLS/' parName '_' runName '_' normType '_harmonics'], 'transitions');
end

%% Plotting
if strcmp(transPlotOrNot, 'yes')
    for iCase = 1:size(transitions,2) % for each type of transition
        
        if isempty(transitions(iCase).f1)
            disp(['no instances of' caseType(iCase)])
            
        elseif size(nonzeros(caseTrials(:,iCase)),1) < minInstances
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

