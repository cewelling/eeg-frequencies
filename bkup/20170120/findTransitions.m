function [ lowToHigh, highToLow, l2hTrials, h2lTrials, l2hMids, h2lMids ] = findTransitions( epochs )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% load parameters
analysisParams

%Find and store transitions between percepts that meet criteria
    
lowToHigh = []; l2hTrials = []; l2hMids = [];
highToLow = []; h2lTrials = []; h2lMids = [];
percType = 'mixed';
epoch1end = NaN;
gap = 0;

for iEpoch = 1:(size(epochs,1))
    thisEpoch = epochs(iEpoch, :);
    if iEpoch > 1
        prevEpoch = epochs(iEpoch - 1, :);
        gap = gap + (thisEpoch(2) - prevEpoch(3));
    end
    epochLength = thisEpoch(3) - thisEpoch(2);
    switch percType
        %------------------------------------------------------------------
        case 'mixed'
            gap = 0;
            if thisEpoch(5) == -1
                if epochLength > domMin
                    percType = 'low';
                    epoch1end = thisEpoch(3); % this epoch marked as pre-transition epoch
                else
                    percType = 'mixed';
                end
            elseif thisEpoch(5) == 1
                if epochLength > domMin
                    percType = 'high';
                    epoch1end = thisEpoch(3); % this epoch marked as pre-transition epoch
                else
                    percType = 'mixed'; % epoch not long enough, pre-transition epoch still unknown
                end
            elseif thisEpoch(5) == 0
                percType = 'mixed';
            end
            %------------------------------------------------------------------
        case 'high'
            if gap + epochLength < gapMax && gap > 0
                percType = 'high';
                gap = gap + epochLength;
            elseif thisEpoch(5) == -1
                if epochLength > domMin
                    if gap < gapMax && thisEpoch(2) - epoch1end > 0
                        transPt = mean([thisEpoch(2) epoch1end]);
                        highToLow = [highToLow; transPt];
                        h2lTrials = [h2lTrials; thisEpoch(1)];
                        mid = thisEpoch(2) - epoch1end;
                        h2lMids = [h2lMids; mid];
                    end
                    percType = 'low';
                    epoch1end = thisEpoch(3);
                else
                    percType = 'mixed';
                end
                gap = 0;
            elseif thisEpoch(5) == 1
                if epochLength > domMin
                    percType = 'high';
                    epoch1end = thisEpoch(3);
                else
                    percType = 'mixed';
                end
                gap = 0;
            elseif thisEpoch(5) == 0
                if epochLength > mixedMin
                    if gap < gapMax && thisEpoch(2) - epoch1end > 0
                        percType = 'mixedFromHigh';
                    else
                        percType = 'mixed';
                    end
                else
                    percType = 'mixed';
                end
                gap = 0;
            end
            %------------------------------------------------------------------
        case 'low'
            if gap + epochLength < gapMax && gap > 0
                percType = 'low';
                gap = gap + epochLength;
            elseif thisEpoch(5) == -1
                if epochLength > domMin
                    percType = 'low';
                    epoch1end = thisEpoch(3);
                else
                    percType = 'mixed';
                end
                gap = 0;
            elseif thisEpoch(5) == 1
                if epochLength > domMin
                    if gap < gapMax && thisEpoch(2) - epoch1end > 0
                        transPt = mean([thisEpoch(2) epoch1end]);
                        lowToHigh = [lowToHigh; transPt];
                        l2hTrials = [l2hTrials; thisEpoch(1)];
                        mid = thisEpoch(2) - epoch1end;
                        l2hMids = [l2hMids; mid];
                    end
                    percType = 'high';
                    epoch1end = thisEpoch(3);
                else
                    percType = 'mixed';
                end
                gap = 0;
            elseif thisEpoch(5) == 0
                if epochLength > mixedMin
                    if gap < gapMax && thisEpoch(2) - epoch1end > 0
                        percType = 'mixedFromLow';
                    else
                        percType = 'mixed';
                    end
                else
                    percType = 'mixed';
                end
                gap = 0;
            end
            %------------------------------------------------------------------
        case 'mixedFromHigh'
            if gap + epochLength < gapMax && gap > 0
                percType = 'mixedFromHigh';
                gap = gap + epochLength;
            elseif thisEpoch(5) == -1
                if epochLength > domMin
                    if gap < gapMax && thisEpoch(2) - epoch1end > 0
                        transPt = mean([thisEpoch(2) epoch1end]);
                        highToLow = [highToLow; transPt];
                        h2lTrials = [h2lTrials; thisEpoch(1)];
                        mid = thisEpoch(2) - epoch1end;
                        h2lMids = [h2lMids; mid];
                    end
                    percType = 'low';
                    epoch1end = thisEpoch(3);
                else
                    percType = 'mixed';
                end
                gap = 0;
            elseif thisEpoch(5) == 1
                if epochLength > domMin
                    percType = 'high';
                    epoch1end = thisEpoch(3);
                else
                    percType = 'mixed';
                end
                gap = 0;
            elseif thisEpoch(5) == 0
                if epochLength > mixedMin
                    if gap < gapMax && thisEpoch(2) - epoch1end > 0
                        percType = 'mixedFromHigh';
                    else
                        percType = 'mixed';
                    end
                else
                    percType = 'mixed';
                end
                gap = 0;
            end
            %------------------------------------------------------------------
        case 'mixedFromLow'
            if 0 < gap + epochLength < gapMax && gap > 0
                percType = 'mixedFromLow';
                gap = gap + epochLength;
            elseif thisEpoch(5) == -1
                if epochLength > domMin
                    percType = 'low';
                    epoch1end = thisEpoch(3);
                else
                    percType = 'mixed';
                end
                gap = 0;
            elseif thisEpoch(5) == 1
                if epochLength > domMin
                    if gap < gapMax && thisEpoch(2) - epoch1end > 0
                        transPt = mean([thisEpoch(2) epoch1end]);
                        lowToHigh = [lowToHigh; transPt];
                        l2hTrials = [l2hTrials; thisEpoch(1)];
                        mid = thisEpoch(2) - epoch1end;
                        l2hMids = [l2hMids; mid];
                    end
                    percType = 'high';
                    epoch1end = thisEpoch(3);
                else
                    percType = 'mixed';
                end
                gap = 0;
            elseif thisEpoch(5) == 0
                if epochLength > mixedMin
                    if gap < gapMax && thisEpoch(2) - epoch1end > 0
                        percType = 'mixedFromLow';
                    else
                        percType = 'mixed';
                    end
                else
                    percType = 'mixed';
                end
                gap = 0;
            end
    end
end
end

