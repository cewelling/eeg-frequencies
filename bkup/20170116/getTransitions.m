function [ tCount ] = getTransitions( parName, runName, date )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Set-up

% load parameters
analysisParams

% get list of participant's perceptual epochs
[epochs, ~] = getEpochs(parName, runName, date);

% codes
high = 1;
mixed = 0;
low = -1;

tCount = 0;
currentT = 'none';
for iEpoch = 1:(size(epochs,1) - 1) %don't check the transition out of the last epoch
    thisEpoch = epochs(iEpoch, :);
    switch currentT
        case 'fromHigh'
            if thisEpoch(5) == high;
                currentT = 'fromHigh';
            elseif thisEpoch(5) == low;
                currentT = 'fromLow';
                tCount = tCount + 1;
            elseif thisEpoch(5) == mixed;
                currentT = 'mfromHigh';
            end
        case 'fromLow'
            if thisEpoch(5) == high;
                currentT = 'fromHigh';
                tCount = tCount + 1;
            elseif thisEpoch(5) == low;
                currentT = 'fromLow';
            elseif thisEpoch(5) == mixed;
                currentT = 'mfromLow';
            end
        case 'mfromHigh'
            if thisEpoch(5) == high;
                currentT = 'fromHigh';
                tCount = tCount + 1;
            elseif thisEpoch(5) == low;
                currentT = 'fromLow';
                tCount = tCount + 1;
            elseif thisEpoch(5) == mixed;
                currentT = 'mfromHigh';
            end
        case 'mfromLow'
            if thisEpoch(5) == high;
                currentT = 'fromHigh';
                tCount = tCount + 1;
            elseif thisEpoch(5) == low;
                currentT = 'fromLow';
                tCount = tCount + 1;
            elseif thisEpoch(5) == mixed;
                currentT = 'mfromLow';
            end
        case 'none' % have not yet entered first dominant percept
            if thisEpoch(5) == high;
                currentT = 'fromHigh';
            elseif thisEpoch(5) == low;
                currentT = 'fromLow';
            elseif thisEpoch(5) == mixed;
                currentT = 'none';
            end
    end
end

tRate = tCount/trialDur; % transitions per second
save ([indicesDir 'tRate/' parName '_' runName], 'tRate');




% % Transition codes to look for in list of epochs, and their labels
% caseCode = [-1 1; 1 -1; 0 1; 0 -1; 1 0; -1 0];
% caseType = {'low to high'; 'high to low'; 'mixed to high'; 'mixed to low'; 'high to mixed'; 'low to mixed'};
% 
% % Preallocate space
% caseTimes = zeros(50,6);
% caseTrials = zeros(50,6);
% 
% % % Find and store transitions between percepts that meet criteria
% 
% lowHigh = []; lhTrials = [];
% highLow = []; hlTrials = [];
% lowMhigh = []; lmhTrials = [];
% highMlow = []; hmlTrials = [];
% lowMlow = []; lmlTrials = [];
% highMhigh = []; hmhTrials = [];
% % lowToHigh = []; l2hTrials = [];
% % highToLow = []; h2lTrials = [];
% % mixedToHigh = []; m2hTrials = [];
% % mixedToLow = []; m2lTrials = [];
% % highToMixed = []; h2mTrials = [];
% % lowToMixed = []; l2mTrials = [];
% 
% for iEpoch = 1:(size(epochs,1) - 1) %don't check the transition out of the last epoch
%     thisEpoch = epochs(iEpoch, :);
%     nextEpoch = epochs(iEpoch + 1, :);
%     nNextEpoch = epochs(iEpoch + 2, :);
%     if nextEpoch(2) - thisEpoch(3) <= gapMax && nextEpoch(1) == thisEpoch(1) % epochs must be in the same trial
%         if thisEpoch(5) == -1
%             if thisEpoch(3) - thisEpoch(2) >= domMin
%                 if nextEpoch(5) == 1
%                     if nextEpoch(3) - nextEpoch(2) >= domMin
%                         transPt = mean([thisEpoch(3) nextEpoch(2)]);
%                         lowtoHigh = [lowToHigh; transPt];
%                         l2hTrials = [l2hTrials; thisEpoch(1)];
%                     end
%                 elseif nextEpoch(5) == 0
%                     if nextEpoch(3) - nextEpoch(2) >= mixedMin
%                         transPt = mean([thisEpoch(3) nextEpoch(2)]);
%                         lowToMixed = [lowToMixed; transPt];
%                         l2mTrials = [l2mTrials; thisEpoch(1)];
%                     end
%                 end
%             end
%             
%         elseif thisEpoch(5) == 1
%             if thisEpoch(3) - thisEpoch(2) >= domMin
%                 if nextEpoch(5) == -1
%                     if nextEpoch(3) - nextEpoch(2) >= domMin
%                         transPt = mean([thisEpoch(3) nextEpoch(2)]);
%                         highToLow = [highToLow; transPt];
%                         h2lTrials = [h2lTrials; thisEpoch(1)];
%                     end
%                 elseif nextEpoch(5) == 0
%                     if nextEpoch(3) - nextEpoch(2) >= mixedMin
%                         transPt = mean([thisEpoch(3) nextEpoch(2)]);
%                         highToMixed = [highToMixed; transPt];
%                         h2mTrials = [h2mTrials; thisEpoch(1)];
%                     end
%                 end
%             end
%             
%         elseif thisEpoch(5) == 0
%             if thisEpoch(3) - thisEpoch(2) >= mixedMin
%                 if nextEpoch(5) == -1
%                     if nextEpoch(3) - nextEpoch(2) >= domMin
%                         transPt = mean([thisEpoch(3) nextEpoch(2)]);
%                         mixedToLow = [mixedToLow; transPt];
%                         m2lTrials = [m2lTrials; thisEpoch(1)];
%                     end
%                 elseif nextEpoch(5) == 1
%                     if nextEpoch(3) - nextEpoch(2) >= domMin
%                         transPt = mean([thisEpoch(3) nextEpoch(2)]);
%                         mixedToHigh = [mixedToHigh; transPt];
%                         m2hTrials = [m2hTrials; thisEpoch(1)];
%                     end
%                 end
%             end
%         end
%     end
% end
% 
end

