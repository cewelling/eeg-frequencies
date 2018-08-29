function [ crossIndex, dist ] = findCross( f1, f2, startIndex, direction)
%findCross.m: For the average of each list of transitions in tList, finds 
% the nearest index to startIndex where the time series for the two stimulation 
% frequencies (f1 and f2) cross. Searches in the direction specified.

% load parameters
groupAnalysisParams
fullLength = 2*transHalf*sampRate + 1;

% % average transitions for each stim frequency
% f1avg = nanmean(tList.f1, 1);
% f2avg = nanmean(tList.f2, 1);

% establish parsing direction
if strcmp(direction, 'right')
    step = 1;
elseif strcmp(direction, 'left')
    step = -1;
end

% parse transition
topFreq = 'none'; % frequency with the higher amplitude
i = startIndex;
while i >= 1 && i <= length(f1) && i <= length(f2)
    switch topFreq
        case 'none'
            if f1(i) > f2(i)
                topFreq = 'lowFreq';
            elseif f2(i) > f1(i)
                topFreq = 'highFreq';
            else
                topFreq = 'none';
            end
        case 'lowFreq'
            if f1(i) > f2(i)
                topFreq = 'lowFreq';
            elseif f2(i) >= f1(i)
                crossIndex = i;
                dist = abs(startIndex - i);
                return;
            end
        case 'highFreq'
            if f2(i) > f1(i)
                topFreq = 'highFreq';
            elseif f1(i) >= f2(i)
                crossIndex = i;
                dist = abs(startIndex - i);
                return;
            end
    end
    i = i + step;
end

% reached edge of data before detecting a cross
crossIndex = NaN;
dist = NaN;
%error('Cross not detected. Increase transHalf.');

end

