function [ peakIndex ] = findPeak(epoch, rls_data, rls_time)
% [peakIndex] = findPeak(epoch, rls_data, rls_time finds and returns the 
% index of the peak of the designated epoch. First finds the absolute 
% maximum of the epoch and keeps it if it is also a local maximum. If not, 
% climbs in the uphill direction until finding a local maximum. 
% 
% Called from: aggPeaks.m
% Dependencies: analysisParams.m


% load parameters
analysisParams

% get trial and dominant frequency
trial = epoch(1);
domFreq = epoch(5);

% get RLS data from appropriate trial
if domFreq == -1 % low freq dominant
    RLS = rls_data(1).amp{trial};
elseif domFreq == 1 % high freq dominant
    try
    RLS = rls_data(2).amp{trial};
    catch
        disp('hi')
    end
end

% find the index of the beginning and end of the epoch
[~, startPtIndex] = min(abs(rls_time - epoch(2)));
[~, endPtIndex] = min(abs(rls_time - epoch(3)));

%% Alternate Peak picking (just find the center of the reported epoch)
if altPeak == 1
    peakIndex = mean([startPtIndex endPtIndex]);
    return;
end

%% Normal peak picking
% Find the index of the max of the epoch
[~, maxInd] = max(RLS(startPtIndex:endPtIndex));
maxInd = maxInd - 1 + startPtIndex; % index relative to entire time series, not just this epoch 
% figure
% plot(startPtIndex:endPtIndex,RLS(startPtIndex:endPtIndex))
% hold on
% plot(maxInd,RLS(maxInd),'x')
% drawnow

% try maxInd < (length(RLS) - (cutoffWindow / 1000 * 512)) && maxInd > 1;
% catch
%     disp('hi')
% end
% only look for peak if reported epoch max is not on edge of time series
if maxInd < (length(RLS) - (cutoffWindow / 1000 * 512)) && maxInd > 1
    
    % if epoch max is not a local maximum, step up until a local max is found
    if RLS(maxInd) <= RLS(maxInd + 1) % step forward
        i = maxInd;
        step = 1;
        platCount = 0; % counts number points in plateau if top of peak is flat
        found = 0;
        while i < length(RLS) - (cutoffWindow / 1000 * 512)
            [platCount, found] = stepUp(RLS, i, platCount, step, found);
            if found == 1
                peakIndex = i - round(platCount/2);
                break;
            end
            i = i+1;
        end
        
        % no peak found
        if found == 0
            peakIndex = NaN;
        end
        
    elseif RLS(maxInd) < RLS(maxInd - 1) % step backward
        i = maxInd;
        step = -1;
        platCount = 0; % counts number points in plateau if top of peak is flat
        found = 0;
        while i > 1
            [platCount, found] = stepUp(RLS, i, platCount, step, found);
            if found == 1
                peakIndex = i + round(platCount/2);
                break;
            end
            i = i-1;
        end
        
        % no peak found
        if found == 0
            peakIndex = NaN;
        end
        
        % maximum of reported epoch is a local maximum, identify it as the peak
    else
        peakIndex = maxInd;
    end
    
    % maximum of reported epoch is on edge of time series, do not count as a peak
else
    peakIndex = NaN;
end

% Fit a sinusoidal curve to peak, take maximum of fit as center of peak
% if ~isnan(peakIndex)
%     fIndex = peakIndex - .5*512*peakHalf;
%     lIndex = peakIndex + .5*512*peakHalf;
%     if fIndex < 1
%         fIndex = 1;
%     elseif lIndex > length(RLS) - (cutoffWindow / 1000 * 512)
%         lIndex = length(RLS) - (cutoffWindow / 1000 * 512);
%     end
%     peakX = fIndex:lIndex;
%     peakY = RLS(fIndex:lIndex);
%     
%     % determine initial parameters for minimization
%     yu = max(peakY);
%     yl = min(peakY);
%     yr = (yu-yl);
%     yz = peakY-yu+(yr/2);
%     zx = peakX(yz .* circshift(yz,[0,1]) <= 0);
%     per = 2000; %2*mean(diff(zx)); % estimated period
%     ym = mean(peakY);
%     
%     fit = @(b,x) b(1).*(sin(2*pi*x./b(2) + 2*pi/b(3))) + b(4); % function to fit
%     fcn = @(b) sum((fit(b,peakX)  - peakY).^2); % least squares cost function
%     options = optimset('MaxFunEvals', 1e10, 'MaxIter', 1e10);
%     s = fminsearch(fcn, [yr; per; -1; ym], options);
%     
%     figure
%     plot(peakX, peakY, 'b',  peakX, fit(s,peakX), 'r')
%     
%     % find index of maximum of fit
%     [~,forMax] = max(fit(s, peakIndex:lIndex));
%     [~,backMax] = max(fit(s, fIndex:peakIndex));
%     forDiff = abs(forMax - peakIndex);
%     backDiff = abs(backMax - peakIndex);
%     if forDiff < backDiff
%         peakIndex = peakIndex + forMax - 1;
%     else
%         peakIndex = fIndex + backMax - 1;
%     end
% end

end

