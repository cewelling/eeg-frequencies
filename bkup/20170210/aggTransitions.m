function [ transitions ] = aggTransitions(parName, runName, date, EEGfile)
% aggTransitions.m: Find and save all chunks of RLS data corresponding to
% participant-reported transitions
%   Detailed explanation goes here

%% Set-up

% load parameters
analysisParams

% Transition codes to look for in list of epochs, and their labels
caseCode = [-1 1; 1 -1; 0 1; 0 -1; 1 0; -1 0];
caseType = {'low to high'; 'high to low'; 'low to low'; 'high to high'};

% Preallocate space
caseTimes = zeros(50,4);
caseTrials = zeros(50,4);

%% Find all transitions

if ~exist(['tLists/' parName '_' runName '.mat'], 'file')
    % get list of participant's perceptual epochs
    [epochs, ~] = getEpochs(parName, runName, date);
    
    % find transitions based on perceptual epochs
    [l2h, h2l, l2l, h2h] = findTransitions(epochs, parName, runName);
else
    load(['tLists/' parName '_' runName '.mat'])
end

%%% In l2h, h2l, etc.
% 1: transition points
% 2: transition trials
% 3: length of time between end of first and beginning of second dom state
%%% Criteria %%%
% 4: length of shortest dom state
% 5: length of shortest mix state
% 6: length of largest gap

caseTimes(1:length(l2h{1}),1) = l2h{1};
caseTrials(1:length(l2h{2}),1) = l2h{2};
caseCrit{1} = [l2h{4} l2h{5} l2h{6}];

caseTimes(1:length(h2l{1}),2) = h2l{1};
caseTrials(1:length(h2l{2}),2) = h2l{2};
caseCrit{2} = [h2l{4} h2l{5} h2l{6}];

caseTimes(1:length(l2l{1}),3) = l2l{1};
caseTrials(1:length(l2l{2}),3) = l2l{2};
caseCrit{3} = [l2l{4} l2l{5} l2l{6}];

caseTimes(1:length(h2h{1}),4) = h2h{1};
caseTrials(1:length(h2h{2}),4) = h2h{2};
caseCrit{4} = [h2h{4} h2h{5} h2h{6}];

%% Get the RLS time-amplitude data for the times corresponding to the transitions

if ~exist(['rls_data/' parName '_' runName '.mat'], 'file')
    [rls_data, rls_time] = runRLS(parName, runName, date, EEGfile);
else
    load(['rls_data/' parName '_' runName '.mat'])
end

%Normalize the RLS data

normRLS_data = [];
% Get baseline using average amplitude when no stimulus was on screen during a sim trial
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

halfLength = 512*transHalf; % Remember we have 512 data points per second
fullLength = halfLength*2+1; % twice the half length plus the transition time

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
        
        % epoch is cut off by beginning, pad beginning with NaNs
        if (tTimeIndex - halfLength) < 1
            epochRLSf1 = RLSf1(1:tTimeIndex + halfLength);
            epochRLSf2 = RLSf2(1:tTimeIndex + halfLength);
            padSize = fullLength - length(epochRLSf1);
            
            epochRLSf1 = [nan(1,padSize), epochRLSf1];
            epochRLSf2 = [nan(1,padSize), epochRLSf2];
            
            % epoch is cut off by end, pad end with NaNs
        elseif (tTimeIndex + halfLength) > length(RLSf1) - (cutoffWindow / 1000 * 512); % remember end is padded due to shift forward
            epochRLSf1 = RLSf1(tTimeIndex - halfLength : length(RLSf1) - (cutoffWindow / 1000 * 512));
            epochRLSf2 = RLSf2(tTimeIndex - halfLength : length(RLSf2) - (cutoffWindow / 1000 * 512));
            padSize = fullLength - length(epochRLSf1);
            
            epochRLSf1 = [epochRLSf1,nan(1,padSize)];
            epochRLSf2 = [epochRLSf2,nan(1,padSize)];
            
            % epoch is not cut off
        else
            epochRLSf1 = RLSf1(tTimeIndex-halfLength:tTimeIndex + halfLength);
            epochRLSf2 = RLSf2(tTimeIndex-halfLength:tTimeIndex + halfLength);
        end
        
        f1Epochs = [f1Epochs; epochRLSf1];
        f2Epochs = [f2Epochs; epochRLSf2];
        
    end
    
    transitions(iCase).type = caseType(iCase);
    transitions(iCase).f1 = f1Epochs;
    transitions(iCase).f2 = f2Epochs;
    transitions(iCase).crit = caseCrit{iCase};
    
    save(['tListsRLS/' parName '_' runName '_' normType], 'transitions');
    
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

