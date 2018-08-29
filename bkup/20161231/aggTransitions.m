function [ transitions ] = aggTransitions(parName, runName, date, EEGfile)
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here

%% Set-up

% load parameters
analysisParams

% get list of participant's perceptual epochs
[epochs, ~] = getEpochs(parName, runName, date);

% Transition codes to look for in list of epochs, and their labels
caseCode = [-1 1; 1 -1; 0 1; 0 -1; 1 0; -1 0];
caseType = {'low to high'; 'high to low'; 'mixed to high'; 'mixed to low'; 'high to mixed'; 'low to mixed'};

% Preallocate space
caseTimes = zeros(50,6);
caseTrials = zeros(50,6);

%% Find and store transition points according to keypress data

for iCase = 1:size(caseCode,1)
    
    % for the transition case, set the minimum durations
    
    if caseCode(iCase,1) ~= 0
        firstStateMin = domMin;
    else
        firstStateMin = mixedMin;
    end
    
    if caseCode(iCase,2) ~= 0
        secondStateMin = domMin;
    else
        secondStateMin = mixedMin;
    end
    
    caseEpochs = []; %preallocate for accruing center points
    caseEpochTrials = [];
    
    for iEpoch = 1:(size(epochs,1) - 1) %don't check the transition out of the last epoch

        % for each line in epochs, check to see if epoch meets criteria:
        % column 5 contains "reporting fast, slow, or mixed percept"
        if epochs(iEpoch, 5) == caseCode(iCase,1) && ...  %the first of the two states is correct
                epochs(iEpoch, 3) - epochs(iEpoch, 2) > firstStateMin && ...  %the first of the two states is sufficiently long
                epochs(iEpoch+1, 5) == caseCode(iCase,2) && ...  %the second state is correct
                epochs(iEpoch+1, 3) - epochs(iEpoch+1, 2) > secondStateMin && ... %the second state is sufficiently long
                0 < (epochs(iEpoch+1,2) - epochs(iEpoch,3)) && ...
                (epochs(iEpoch+1,2) - epochs(iEpoch,3)) < gapMax %the time between reports is sufficiently small, and not negative
            
            centerPt = mean([epochs(iEpoch,3) epochs(iEpoch+1,2)]); % if the epoch transition meets criteria, store the center time point
            caseEpochs = [caseEpochs; centerPt];
            caseEpochTrials = [caseEpochTrials; epochs(iEpoch,1)];
            
        end
        
    end
    
    caseTimes(1:length(caseEpochs),iCase) = caseEpochs;
    caseTrials(1:length(caseEpochs),iCase) = caseEpochTrials;
    
end

%% Get the RLS time-amplitude data for the times corresponding to the transitions

[rls_data rls_time] = runRLS(parName, runName, date, EEGfile);

%Normalize the RLS data

normRLS_data = [];

for iTrial = 1:size(rls_data(1).amp, 2)
    
    for iFreq = 1:length(rls_data);
        
        RLSamp = (rls_data(iFreq).amp{iTrial})';
        
        % take normalizing/standardizing step
        if strcmp(normType, 'mean')
            normRLSamp = RLSamp - mean(RLSamp);
        elseif strcmp(normType, 'norm')
            normRLSamp = 2*(RLSamp - min(RLSamp))/(max(RLSamp) - min(RLSamp)) - 1;
        elseif strcmp(normType, 'z')
            normRLSamp = zscore(RLSamp);
        elseif strcmp(normType, 'none')
            normRLSamp = RLSamp;
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
        elseif (tTimeIndex + halfLength) > length(RLSf1) - (timeWindow / 1000 * 512); % remember end is padded due to shift forward
            epochRLSf1 = RLSf1(tTimeIndex - halfLength : length(RLSf1) - (timeWindow / 1000 * 512));
            epochRLSf2 = RLSf2(tTimeIndex - halfLength : length(RLSf2) - (timeWindow / 1000 * 512));
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

