function [ epochEntry ] = handleEpoch(keyPressData, eLastRow, startTime)
% [epochEntry] = handleEpoch(keyPressData, eLastRow, startTime) returns an 
% epochEntry (a vector) with the following elements: 
% 1. the trial in which the epoch took place
% 2. start time of the epoch
% 3. button pressed (100 for left, 102 for right, 104 for up)
% 4. dominant frequency (1 for high, -1 for low, 0 for mixed)
%
% Arguments:
% keyPressData is the matrix of behavioral data produced by the runScript
% eLastRow is the last row of the epoch we are storing. 
% startTime is the start time of the epoch
%
% eLastRow and startTime are determined in getEpochs.m
%
% Called from: getEpochs.m
% Requires: analysisParams.m
% Created by Alina Spiegel 11-18-2016
%--------------------------------------------------------------------------

%load params
analysisParams

endTime = keyPressData(eLastRow - 1, 3);
currTrial = keyPressData(eLastRow, 1);
eLastRow = eLastRow - 1; % last row of the epoch we are storing
if keyPressData(eLastRow, 4) ~= keyPressData(eLastRow, 5) && keyPressData(eLastRow, 2) == leftArrow
    domFreq = 1; %high
elseif keyPressData(eLastRow, 4) ~= keyPressData(eLastRow, 5) && keyPressData(eLastRow, 2) == rightArrow
    domFreq = -1; %low
elseif keyPressData(eLastRow, 4) == keyPressData(eLastRow, 5) && keyPressData(eLastRow, 2) == leftArrow
    domFreq = -1;
elseif keyPressData(eLastRow, 4) == keyPressData(eLastRow, 5) && keyPressData(eLastRow, 2) == rightArrow
    domFreq = 1;
else % mixed epoch
    domFreq = 0;
end

epochEntry = [currTrial startTime endTime keyPressData(eLastRow, 2) domFreq];

end

