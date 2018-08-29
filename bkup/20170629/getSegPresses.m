function [ segPressMat ] = getSegPresses( parName, runName, tTrialVect, tPointVect, rivLat)
%getSegPresses Summary of this function goes here
%   Outputs a matrix of button presses corresponding to the designated transitions
% Each row is a transition, each column is a time point within the
% transition (transition happens at time 0)

% load parameters
analysisParams

if strfind(runName, 'Sim')
    % load sim schedules
    load(['simSchedules/' parName '_' runName '.mat'])
    keyPressRef = simSchedule;
    latency = 0;
else
    % load epochs
    load(['epochs/' parName '_' runName '.mat'])
    % remember that tPoints have been shifted by 1 second, but epochs haven't **
    keyPressRef = epochs;
    latency = rivLat; %0.566; % latency for rivalry runs
end

% time matrix for transitions
tRange = -transHalf:1/sampRate:transHalf;
tPoints = repmat(tPointVect,[1,size(tRange, 2)]);
tRanges = repmat(tRange, [size(tPoints, 1), 1]);
tTimeMat = tPoints + tRanges;

% account for this: tPoints have been shifted back by 1 second, but epochs haven't **
tTimeMat = tTimeMat + 1;

% space to store button press information
segPressMat = [];

% iterate through epochs, identify time points that correspond to epochs
for thisTrial = min(keyPressRef(:,1)):max(keyPressRef(:,1))
    trialEpochs = keyPressRef(keyPressRef(:,1) == thisTrial,:);
    trialTimes = tTimeMat(tTrialVect == thisTrial,:);
    
    % space to store button press information
    trialPressMat = nan(size(trialTimes, 1), size(trialTimes, 2));
    
    for iEpoch = 1:size(trialEpochs,1)
        
        % mark trial times with the appropriate button presses
        trialPressMat(trialTimes > (trialEpochs(iEpoch,2) - latency) & trialTimes < (trialEpochs(iEpoch, 3) - latency)) = trialEpochs(iEpoch, 5);
    end
    
    % combine trial matrices into one run matrix
    segPressMat = [segPressMat; trialPressMat];
end

end

