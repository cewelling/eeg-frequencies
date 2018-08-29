function [epochs, simSchedule] = getEpochs( parName, runName, date)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

% load parameters
analysisParams
    
% load key press data file
keyPressData = load([keyPressDir parName '_' runName '_' date '.txt']);

%% Parse key press data

epochs = [];
startTime = NaN;
currEpoch = 'none';
for currRow = 1:size(keyPressData,1)
    switch currEpoch
        case 'none'
            if keyPressData(currRow, 2) == leftArrow
                currEpoch = 'left';
                startTime = keyPressData(currRow, 3);
            elseif keyPressData(currRow, 2) == rightArrow
                currEpoch = 'right';
                startTime = keyPressData(currRow, 3);
            elseif keyPressData(currRow, 2) == upArrow
                currEpoch = 'mixed';
                startTime = keyPressData(currRow, 3);
            elseif keyPressData(currRow, 2) == 0
                currEpoch = 'none';
            end
            
        case 'left'
            if keyPressData(currRow, 2) == leftArrow
                currEpoch = 'left';
            elseif keyPressData(currRow, 2) == rightArrow
                epochEntry = handleEpoch(keyPressData, currRow-1, startTime);
                % store epoch
                epochs = [epochs ; epochEntry];
                % reset start time
                startTime = keyPressData(currRow, 3);
                currEpoch = 'right';
            elseif keyPressData(currRow, 2) == upArrow
                epochEntry = handleEpoch(keyPressData, currRow-1, startTime);
                % store epoch
                epochs = [epochs ; epochEntry];
                % reset start time
                startTime = keyPressData(currRow, 3);
                currEpoch = 'mixed';
            elseif keyPressData(currRow, 2) == 0
                epochEntry = handleEpoch(keyPressData, currRow-1, startTime);
                % store epoch
                epochs = [epochs ; epochEntry];
                currEpoch = 'none';               
            end
                           
        case 'right'
            if keyPressData(currRow, 2) == rightArrow
                currEpoch = 'right';
            elseif keyPressData(currRow, 2) == leftArrow
                epochEntry = handleEpoch(keyPressData, currRow-1, startTime);
                % store epoch
                epochs = [epochs ; epochEntry];
                % reset start time
                startTime = keyPressData(currRow, 3);
                currEpoch = 'left';
            elseif keyPressData(currRow, 2) == upArrow
                epochEntry = handleEpoch(keyPressData, currRow-1, startTime);
                % store epoch
                epochs = [epochs ; epochEntry];
                % reset start time
                startTime = keyPressData(currRow, 3);
                currEpoch = 'mixed';
            elseif keyPressData(currRow, 2) == 0
                epochEntry = handleEpoch(keyPressData, currRow-1, startTime);
                % store epoch
                epochs = [epochs ; epochEntry];
                currEpoch = 'none';
            end
            
        case 'mixed'
            if keyPressData(currRow, 2) == upArrow
                currEpoch = 'mixed';
            elseif keyPressData(currRow, 2) == leftArrow
                epochEntry = handleEpoch(keyPressData, currRow-1, startTime);
                % store epoch
                epochs = [epochs ; epochEntry];
                % reset start time
                startTime = keyPressData(currRow, 3);
                currEpoch = 'left';
            elseif keyPressData(currRow, 2) == rightArrow
                epochEntry = handleEpoch(keyPressData, currRow-1, startTime);
                % store epoch
                epochs = [epochs ; epochEntry];
                % reset start time
                startTime = keyPressData(currRow, 3);
                currEpoch = 'right';
            elseif keyPressData(currRow, 2) == 0
                epochEntry = handleEpoch(keyPressData, currRow-1, startTime);
                % store epoch
                epochs = [epochs ; epochEntry];
                currEpoch = 'none';
            end
    end
end

%% Account for mistakes (lost EEG data)
if strcmp(parName, 'cumulus05') && strcmp(runName, 'marzRival1')
    % account for loss of first two trials of this run in EEG data
    lastIndex = find(epochs(:,1) < 3, 1, 'last');
    epochs = epochs(lastIndex + 1:end, :);
end

if strcmp(parName, 'stratus135') && strcmp(runName, 'dartSim1')
    % account for loss of first trial of this run in EEG data
    lastIndex = find(epochs(:,1) < 2, 1, 'last');
    epochs = epochs(lastIndex + 1:end, :);
end

if strcmp(parName, 'cumulus01') && strcmp(runName, 'dartSim3')
    % account for loss of first trial of this run in EEG data
    lastIndex = find(epochs(:,1) < 2, 1, 'last');
    epochs = epochs(lastIndex + 1:end, :);
end

if strcmp(parName, 'cumulus01') && strcmp(runName, 'marzRival3')
    % account for loss of first two trials of this run in EEG data
    lastIndex = find(epochs(:,1) < 3, 1, 'last');
    epochs = epochs(lastIndex + 1:end, :);
end

save(['epochs/' parName '_' runName], 'epochs');

%% Parse simulation schedule

if size(keyPressData, 2) == 6 % if data has presentation schedule column
    simSchedule = [];
    startTime = NaN;
    currEpoch = 'none';
    for currRow = 1:size(keyPressData,1)
        switch currEpoch
            case 'none'
                if keyPressData(currRow, 6) == leftCue % saw green
                    currEpoch = 'left';
                    startTime = keyPressData(currRow, 3);
                elseif keyPressData(currRow, 6) == rightCue % saw red
                    currEpoch = 'right';
                    startTime = keyPressData(currRow, 3);
                elseif keyPressData(currRow, 6) == 0
                    currEpoch = 'none';
                end
                
            case 'left'
                if keyPressData(currRow, 6) == leftCue
                    currEpoch = 'left';
                elseif keyPressData(currRow, 6) == rightCue
                    epochEntry = handleEpoch(keyPressData, currRow, startTime);
                    % store epoch
                    simSchedule = [simSchedule ; epochEntry];
                    % reset start time
                    startTime = keyPressData(currRow, 3);
                    currEpoch = 'right';
                elseif keyPressData(currRow, 6) == 0
                    epochEntry = handleEpoch(keyPressData, currRow, startTime);
                    % store epoch
                    simSchedule = [simSchedule ; epochEntry];
                    currEpoch = 'none';
                end
                
            case 'right'
                if keyPressData(currRow, 6) == rightCue
                    currEpoch = 'right';
                elseif keyPressData(currRow, 6) == leftCue
                    epochEntry = handleEpoch(keyPressData, currRow, startTime);
                    % store epoch
                    simSchedule = [simSchedule ; epochEntry];
                    % reset start time
                    startTime = keyPressData(currRow, 3);
                    currEpoch = 'left';
                elseif keyPressData(currRow, 6) == 0
                    epochEntry = handleEpoch(keyPressData, currRow, startTime);
                    % store epoch
                    simSchedule = [simSchedule ; epochEntry];
                    currEpoch = 'none';
                end
        end
    end
else
    simSchedule = NaN;
end

end

