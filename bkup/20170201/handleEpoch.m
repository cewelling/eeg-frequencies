function [ epochEntry ] = handleEpoch(keyPressData, eLastRow, startTime)
% Description: Create an epochEntry including the trial in which the epoch 
%              took place, start/end time of the epoch, button pressed, 
%              and dominant frequency.
%
% Called from: getEpochs.m
% Requires: analysisParams.m
% Created by Alina Spiegel 11-18-2016
%--------------------------------------------------------------------------

%load params
analysisParams

endTime = keyPressData(eLastRow - 1, 3);
if endTime - startTime >= minEpochDur
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
else
    epochEntry = [];
end

end

