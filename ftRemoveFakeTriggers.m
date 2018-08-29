function [ newPreProc_data ] = ftRemoveFakeTriggers( preProc_data, EEGfile, numTrials )
% ftRemoveFakeTriggers returns pre-processed data with only "real" trials.
% A mistake in our runScripts was sending some misplaced trial-marking
% signals (triggers) to the EEG computer. This function removes the
% 'trials' created by those misplaced triggers in preProc_data.

thisData = preProc_data;


% stratus138
if ~isempty(strfind(EEGfile, 'stratus138_dartRival1')) || ~isempty(strfind(EEGfile, 'stratus138_dartSim3')) ...
        || ~isempty(strfind(EEGfile, 'stratus138_marzRival1'))
    thisData = removeTrials(thisData, numTrials, 3);
elseif ~isempty(strfind(EEGfile, 'stratus138_dart')) || ~isempty(strfind(EEGfile, 'stratus138_marzRival2'))
    thisData = removeTrials(thisData, numTrials, 2, 4);
elseif ~isempty(strfind(EEGfile, 'stratus138'))
    thisData = removeTrials(thisData, numTrials, 2);
end

% cumulus15
if ~isempty(strfind(EEGfile, 'cumulus15_dartSim1')) || ~isempty(strfind(EEGfile, 'cumulus15_dartSim2')) ...
        || ~isempty(strfind(EEGfile, 'cumulus15_marzRival3'))
    thisData = removeTrials(thisData, numTrials, 3);
elseif ~isempty(strfind(EEGfile, 'cumulus15_dartRival1')) || ~isempty(strfind(EEGfile, 'cumulus15_marzRival'))
    thisData = removeTrials(thisData, numTrials, 2, 4);
elseif ~isempty(strfind(EEGfile, 'cumulus15'))
    thisData = removeTrials(thisData, numTrials, 2);
end

% cumulus05
if ~isempty(strfind(EEGfile, 'cumulus05_dartRival2')) || ~isempty(strfind(EEGfile, 'cumulus05_dartSim2'))
    thisData = removeTrials(thisData, numTrials, 2, 4);
elseif ~isempty(strfind(EEGfile, 'cumulus05_dartSim3')) || ~isempty(strfind(EEGfile, 'cumulus05_marzRival2')) ...
        || ~isempty(strfind(EEGfile, 'cumulus05_marzRival3'))
    thisData = removeTrials(thisData, numTrials, 2);
end

% cumulus06
if ~isempty(strfind(EEGfile, 'cumulus06_marzRival1'))
    thisData = removeTrials(thisData, numTrials, 3);
elseif ~isempty(strfind(EEGfile, 'cumulus06_dartRival1')) || ~isempty(strfind(EEGfile, 'cumulus06_dartSim1')) ...
        || ~isempty(strfind(EEGfile, 'cumulus06_dartSim3')) || ~isempty(strfind(EEGfile, 'cumulus06_marzRival'))
    thisData = removeTrials(thisData, numTrials, 2, 4);
elseif ~isempty(strfind(EEGfile, 'cumulus06'))
    thisData = removeTrials(thisData, numTrials, 2);
end

% cumulus28
if ~isempty(strfind(EEGfile, 'cumulus28_dartRival1'))
    thisData = removeTrials(thisData, numTrials, 1);
end

newPreProc_data = thisData;

end

