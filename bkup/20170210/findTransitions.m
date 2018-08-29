function [ l2h, h2l, l2l, h2h ] = findTransitions( epochs, parName, runName )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% load parameters
analysisParams

%Find and store transitions between percepts that meet criteria
    
lowToHigh = []; l2hTrials = []; l2hMids = []; l2hDomL = []; l2hMixedL = []; l2hGapL = []; 
highToLow = []; h2lTrials = []; h2lMids = []; h2lDomL = []; h2lMixedL = []; h2lGapL = [];
lowToLow = []; l2lTrials = []; l2lMids = []; l2lDomL = []; l2lMixedL = []; l2lGapL = [];
highToHigh = []; h2hTrials = []; h2hMids = []; h2hDomL = []; h2hMixedL = []; h2hGapL = [];
% transPt 
% Trials: trial it came from
% Mids: length of space between dominant states
% DomL: length of shortest dom state
% MixedL: length of mixed state
% GapL: length of biggest gap (participant did not have a key down)

percType = 'mixed';
epoch1end = NaN;
epoch1length = NaN;
preMixedGap = [];
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
                percType = 'low';
                epoch1end = thisEpoch(3); % this epoch marked as pre-transition epoch
                epoch1length = epochLength;
            elseif thisEpoch(5) == 1
                percType = 'high';
                epoch1end = thisEpoch(3); % this epoch marked as pre-transition epoch
                epoch1length = epochLength;
            elseif thisEpoch(5) == 0
                percType = 'mixed';
            end
            %------------------------------------------------------------------
        case 'high'
            if thisEpoch(5) == -1                
                if thisEpoch(2) - epoch1end > 0
                    transPt = mean([thisEpoch(2) epoch1end]);
                    highToLow = [highToLow; transPt];
                    h2lTrials = [h2lTrials; thisEpoch(1)];
                    mid = thisEpoch(2) - epoch1end;
                    h2lMids = [h2lMids; mid];
                    h2lDomL = [h2lDomL; min([epochLength epoch1length])];
                    h2lMixedL = [h2lMixedL; NaN];
                    h2lGapL = [h2lGapL; gap];
                end
                percType = 'low';
                epoch1end = thisEpoch(3);
                epoch1length = epochLength;
                gap = 0;
            elseif thisEpoch(5) == 1
                percType = 'high';
                if gap > gapMax
                    epoch1end = thisEpoch(3);
                    epoch1length = epochLength;
                    gap = 0;
                end
            elseif thisEpoch(5) == 0
                if thisEpoch(2) - epoch1end > 0
                    percType = 'mixedFromHigh';
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
                percType = 'low';
                if gap > gapMax
                    epoch1end = thisEpoch(3);
                    epoch1length = epochLength;
                    gap = 0;
                end
            elseif thisEpoch(5) == 1
                if thisEpoch(2) - epoch1end > 0
                    transPt = mean([thisEpoch(2) epoch1end]);
                    lowToHigh = [lowToHigh; transPt];
                    l2hTrials = [l2hTrials; thisEpoch(1)];
                    mid = thisEpoch(2) - epoch1end;
                    l2hMids = [l2hMids; mid];
                    l2hDomL = [l2hDomL; min([epochLength epoch1length])];
                    l2hMixedL = [l2hMixedL; NaN];
                    l2hGapL = [l2hGapL; gap];
                end
                percType = 'high';
                epoch1end = thisEpoch(3);
                epoch1length = epochLength;
                gap = 0;
            elseif thisEpoch(5) == 0
                if epochLength > mixedMin
                    if thisEpoch(2) - epoch1end > 0
                        percType = 'mixedFromLow';
                        preMixedGap = gap;
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
            if thisEpoch(5) == -1
                if thisEpoch(2) - epoch1end > 0
                    transPt = mean([thisEpoch(2) epoch1end]);
                    highToLow = [highToLow; transPt];
                    h2lTrials = [h2lTrials; thisEpoch(1)];
                    mid = thisEpoch(2) - epoch1end;
                    h2lMids = [h2lMids; mid];
                    h2lDomL = [h2lDomL; min([epochLength epoch1length])];
                    h2lMixedL = [h2lMixedL; prevEpoch(3) - prevEpoch(2)];
                    h2lGapL = [h2lGapL; max([gap preMixedGap])];
                end
                percType = 'low';
                epoch1end = thisEpoch(3);
                preMixedGap = [];
                gap = 0;
            elseif thisEpoch(5) == 1
                if thisEpoch(2) - epoch1end > 0
                    revPt = mean([thisEpoch(2) epoch1end]);
                    highToHigh = [highToHigh; revPt];
                    h2hTrials = [h2hTrials; thisEpoch(1)];
                    mid = thisEpoch(2) - epoch1end;
                    h2hMids = [h2hMids; mid];
                    h2hDomL = [h2hDomL; min([epochLength epoch1length])];
                    h2hMixedL = [h2hMixedL; prevEpoch(3) - prevEpoch(2)];
                    h2hGapL = [h2hGapL; max([gap preMixedGap])];
                end
                percType = 'high';
                epoch1end = thisEpoch(3);
                epoch1length = epochLength;
                preMixedGap = [];
                gap = 0;
            elseif thisEpoch(5) == 0
                if thisEpoch(2) - epoch1end > 0
                    percType = 'mixedFromHigh';
                    preMixedGap = [preMixedGap gap];
                else
                    percType = 'mixed';
                end
                gap = 0;
            end
            %------------------------------------------------------------------
        case 'mixedFromLow'
            if thisEpoch(5) == -1
                if thisEpoch(2) - epoch1end > 0
                    revPt = mean([thisEpoch(2) epoch1end]);
                    lowToLow = [lowToLow; revPt];
                    l2lTrials = [l2lTrials; thisEpoch(1)];
                    mid = thisEpoch(2) - epoch1end;
                    l2lMids = [l2lMids; mid];
                    l2lDomL = [l2lDomL; min([epochLength epoch1length])];
                    l2lMixedL = [l2lMixedL; prevEpoch(3) - prevEpoch(2)];
                    l2lGapL = [l2lGapL; max([gap preMixedGap])];
                end
                percType = 'low';
                epoch1end = thisEpoch(3);
                epoch1length = epochLength;
                preMixedGap = [];
                gap = 0;
            elseif thisEpoch(5) == 1
                if thisEpoch(2) - epoch1end > 0
                    transPt = mean([thisEpoch(2) epoch1end]);
                    lowToHigh = [lowToHigh; transPt];
                    l2hTrials = [l2hTrials; thisEpoch(1)];
                    mid = thisEpoch(2) - epoch1end;
                    l2hMids = [l2hMids; mid];
                    l2hDomL = [l2hDomL; min([epochLength epoch1length])];
                    l2hMixedL = [l2hMixedL; prevEpoch(3) - prevEpoch(2)];
                    l2hGapL = [l2hGapL; max([gap preMixedGap])];   
                end
                percType = 'high';
                epoch1end = thisEpoch(3);
                epoch1length = epochLength;
                preMixedGap = [];
                gap = 0;
            elseif thisEpoch(5) == 0
                if thisEpoch(2) - epoch1end > 0
                    percType = 'mixedFromLow';
                    preMixedGap = [preMixedGap gap];
                else
                    percType = 'mixed';
                end
                gap = 0;
            end
    end
end

l2h = {lowToHigh, l2hTrials, l2hMids, l2hDomL, l2hMixedL, l2hGapL};
h2l = {highToLow, h2lTrials, h2lMids, h2lDomL, h2lMixedL, h2lGapL};
l2l = {lowToLow, l2lTrials, l2lMids, l2lDomL, l2lMixedL, l2lGapL}; 
h2h = {highToHigh, h2hTrials, h2hMids, h2hDomL, h2hMixedL, h2hGapL}; 

save(['tLists/' parName '_' runName], 'l2h', 'h2l', 'l2l', 'h2h') 

end

