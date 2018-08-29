function [ l2h, h2l, l2l, h2h ] = parseSwitchData( parName, runName )
% parseSwitchData(parName, runName) parses Caroline's switch data matrices
% and saves a cell array corresponding to each type of transition.  
% [l2h, h2l, l2l, h2h] = parseSwitchData(parName, runName) returns the 
% aforementioned cell arrays. Each cell corresponds to a transition
% property as follows:
%
% 1. Time point corresponding to the transition 
% 2. Trial containing the transition
% 3. Time from the end of the 1st dominant epoch to beginning of the 2nd
% 4. Length of 1st dominant epoch
% 5. Length of 2nd dominant epoch
% 6. Length of mixed epoch
% 7. Length of largest gap 
%
% *** SHIFTS transitions 1 second back to account for 1 second shift of RLS 
% data!!! ***
%
% Called from: aggTransitions.m
% Dependencies: analysisParams.m

% load parameters
analysisParams

% Create space to store transitions and their properties
lowToHigh = []; l2hTrials = []; l2hMids = []; l2hDom1 = []; l2hDom2 = []; l2hMixedL = []; l2hGapL = [];
highToLow = []; h2lTrials = []; h2lMids = []; h2lDom1 = []; h2lDom2 = []; h2lMixedL = []; h2lGapL = [];
lowToLow = []; l2lTrials = []; l2lMids = []; l2lDom1 = []; l2lDom2 = []; l2lMixedL = []; l2lGapL = [];
highToHigh = []; h2hTrials = []; h2hMids = []; h2hDom1 = []; h2hDom2 = []; h2hMixedL = []; h2hGapL = [];

% get trial numbers for switchData file identification
runNum = str2double(runName(end));
trialNums = (runNum - 1) * numTrials + 1 : runNum * numTrials;

% account for lost trials (due to data collection mistakes)
if strcmp(parName, 'cumulus05') && strcmp(runName, 'marzRival1')
    % account for loss of first two trials of this run in EEG data
    trialNums = 3:6;
end

if strcmp(parName, 'stratus135') && strcmp(runName, 'dartSim1')
    % account for loss of first trial of this run in EEG data
    trialNums = 2:6;
end

if strcmp(parName, 'cumulus01') && strcmp(runName, 'dartSim3')
    % account for loss of first trial of this run in EEG data
    trialNums = 14:18;
end

if strcmp(parName, 'cumulus01') && strcmp(runName, 'marzRival3')
    % account for loss of first two trials of this run in EEG data
    trialNums = 15:18;
end

% parse appropriate trials
for iTrial = trialNums
    
    % load switch data for this trial
    if exist(['../Behavior/switchData/' parName '_' num2str(iTrial) '_' runName(1:end - 1) '_switchData.mat'], 'file')
        load(['../Behavior/switchData/' parName '_' num2str(iTrial) '_' runName(1:end - 1) '_switchData.mat'])
    else
        error(['No switch data for ' parName '_' runName(1:end - 1) ' trial ' num2str(iTrial)])
    end
    
    % populate transition lists using switch data
    for i = 1:size(switchData, 1)
        
        if switchData(i, 10) == 1 % indicates a switch
            
            % Low to High
            if switchData(i, 4) == 10 % transition starts with low frequency
                
                % initialize variables
                mixedLength = 0;
                gapLength = 0;
                
                % get transition properties
                if switchData(i + 1, 4) == 11
                    transRow = i + 1;
                    gapLength = switchData(i + 1, 5) - switchData(i, 6);
                elseif switchData(i + 2, 4) == 11
                    transRow = i + 2;
                    mixedLength = switchData(i + 1, 7); % middle must be mixed
                    gap1 = switchData(i + 1, 5) - switchData(i, 6);
                    gap2 = switchData(i + 2, 5) - switchData(i + 1, 6);
                    gapLength = max([gap1 gap2]);
                end
                
                % store transition
                lowToHigh = [lowToHigh; mean([switchData(transRow, 5) switchData(i, 6)]) - 1]; % subtract 1 to account for 1 second of EEG data cut off at beginning (onset response)
                l2hTrials = [l2hTrials; iTrial - (numTrials * (runNum - 1))];
                l2hMids = [l2hMids; switchData(transRow, 5) - switchData(i, 6)];
                l2hDom2 = [l2hDom2; switchData(transRow, 7)];
                l2hDom1 = [l2hDom1; switchData(i, 7)];
                l2hMixedL = [l2hMixedL; mixedLength];
                l2hGapL = [l2hGapL; gapLength];
                
                % High to Low
            elseif switchData(i, 4) == 11 % transition starts with high frequency
                
                % initialize variables
                mixedLength = 0;
                gapLength = 0;
                
                % get transition properties
                if switchData(i + 1, 4) == 10
                    transRow = i + 1;
                    gapLength = switchData(i + 1, 5) - switchData(i, 6);
                elseif switchData(i + 2, 4) == 10
                    transRow = i + 2;
                    mixedLength = switchData(i + 1, 7); % middle must be mixed
                    gap1 = switchData(i + 1, 5) - switchData(i, 6);
                    gap2 = switchData(i + 2, 5) - switchData(i + 1, 6);
                    gapLength = max([gap1 gap2]);
                end
                
                % store transition
                highToLow = [highToLow; mean([switchData(transRow, 5) switchData(i, 6)]) - 1]; % subtract 1 to account for 1 second of EEG data cut off at beginning (onset response)
                h2lTrials = [h2lTrials; iTrial - (numTrials * (runNum - 1))];
                h2lMids = [h2lMids; switchData(transRow, 5) - switchData(i, 6)];
                h2lDom2 = [h2lDom2; switchData(transRow, 7)];
                h2lDom1 = [h2lDom1; switchData(i, 7)];
                h2lMixedL = [h2lMixedL; mixedLength];
                h2lGapL = [h2lGapL; gapLength];
                
            end
            
        elseif switchData(i, 10) == 2
            
            % Low to low
            if switchData(i, 4) == 10 % transition starts with low frequency
                
                % initialize variables
                mixedLength = 0;
                gapLength = 0;
                
                % get transition properties
                if switchData(i + 1, 4) == 10
                    transRow = i + 1;
                    gapLength = switchData(i + 1, 5) - switchData(i, 6);
                elseif switchData(i + 2, 4) == 10
                    transRow = i + 2;
                    mixedLength = switchData(i + 1, 7); % middle must be mixed
                    gap1 = switchData(i + 1, 5) - switchData(i, 6);
                    gap2 = switchData(i + 2, 5) - switchData(i + 1, 6);
                    gapLength = max([gap1 gap2]);
                end
                
                % store transition
                lowToLow = [lowToLow; mean([switchData(transRow, 5) switchData(i, 6)]) - 1]; % subtract 1 to account for 1 second of EEG data cut off at beginning (onset response)
                l2lTrials = [l2lTrials; iTrial - (numTrials * (runNum - 1))];
                l2lMids = [l2lMids; switchData(transRow, 5) - switchData(i, 6)];
                l2lDom2 = [l2lDom2; switchData(transRow, 7)];
                l2lDom1 = [l2lDom1; switchData(i, 7)];
                l2lMixedL = [l2lMixedL; mixedLength];
                l2lGapL = [l2lGapL; gapLength];
                
                % High to high
            elseif switchData(i, 4) == 11 % transition starts with high frequency
                
                % initialize variables
                mixedLength = 0;
                gapLength = 0;
                
                % get transition properties
                if switchData(i + 1, 4) == 11
                    transRow = i + 1;
                    gapLength = switchData(i + 1, 5) - switchData(i, 6);
                elseif switchData(i + 2, 4) == 11
                    transRow = i + 2;
                    mixedLength = switchData(i + 1, 7); % middle must be mixed
                    gap1 = switchData(i + 1, 5) - switchData(i, 6);
                    gap2 = switchData(i + 2, 5) - switchData(i + 1, 6);
                    gapLength = max([gap1 gap2]);
                end
                
                % store transition
                highToHigh = [highToHigh; mean([switchData(transRow, 5) switchData(i, 6)]) - 1]; % subtract 1 to account for 1 second of EEG data cut off at beginning (onset response)
                h2hTrials = [h2hTrials; iTrial - (numTrials * (runNum - 1))];
                h2hMids = [h2hMids; switchData(transRow, 5) - switchData(i, 6)];
                h2hDom2 = [h2hDom2; switchData(transRow, 7)];
                h2hDom1 = [h2hDom1; switchData(i, 7)];
                h2hMixedL = [h2hMixedL; mixedLength];
                h2hGapL = [h2hGapL; gapLength];
                
                
            end
        end
    end
end

l2h = {lowToHigh, l2hTrials, l2hMids, l2hDom1, l2hDom2, l2hMixedL, l2hGapL};
h2l = {highToLow, h2lTrials, h2lMids, h2lDom1, h2lDom2, h2lMixedL, h2lGapL};
l2l = {lowToLow, l2lTrials, l2lMids, l2lDom1, l2lDom2, l2lMixedL, l2lGapL};
h2h = {highToHigh, h2hTrials, h2hMids, h2hDom1, h2hDom2, h2hMixedL, h2hGapL};

save(['transitions/tLists_CERW/' parName '_' runName], 'l2h', 'h2l', 'l2l', 'h2h')

end

