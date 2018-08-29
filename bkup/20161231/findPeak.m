function [ peakIndex ] = findPeak(epoch, rls_data, rls_time)
% Find index of the peak of the designated epoch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Start at reported end of the epoch, crawl back until finding a local maximum greater
%   than or equal to the absolute max of the reported percept 
%   Called from: findPeak.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get trial and dominant frequency
trial = epoch(1);
domFreq = epoch(5);

% get RLS data from appropriate trial
if domFreq == -1 % low freq dominant
    RLS = rls_data(1).amp{trial};
elseif domFreq == 1; % high freq dominant
    RLS = rls_data(2).amp{trial};
end

% find the index of the beginning and end of the epoch
[~, startPtIndex] = min(abs(rls_time - epoch(2)));
[~, endPtIndex] = min(abs(rls_time - epoch(3)));

% Find the max of the epoch
epochMax = max(RLS(startPtIndex:endPtIndex));
% if epochMax == 2.8638
%     disp('hi')
% end

% crawl forward from startpoint, backward from endpoint, look for max, make 
% sure that it is at least as big as the max of the epoch

% starting conditions
i = startPtIndex; j = endPtIndex;
frontStep = 1; backStep = -1;
platCount = 0; % counts number points in plateau if top of peak is flat
forSlope = findStartSlope(RLS, i, frontStep); 
backSlope = findStartSlope(RLS, j, backStep); 

% alternate stepping forward and backward
peakIndex = stepForward(RLS, epochMax, i, j, forSlope, backSlope, platCount, frontStep); 


% if endPtIndex == 3628
%     figure
%     hold on
% end
% platCount = 0; % counts number points in plateau if top of peak is flat
% slope = findStartSlope(RLS, i);
% while(1)   
% %     if i == 1
% %         disp('hi')
%     end
%     switch slope
%         case 'incline'
%             if RLS(i) >= epochMax
%                 if RLS(i) > RLS(i-1)
%                     peakIndex = i + round(platCount/2); % choose middle of plateau
%                     return;
%                 elseif RLS(i) == RLS(i-1)
%                     platCount = platCount + 1; 
%                     slope = 'incline';
%                 end
%             elseif RLS(i) > RLS(i-1)
%                 platCount = 0;
%                 slope = 'decline';
%             else
%                 platCount = 0;
%                 slope = 'incline';
%             end
%                     
%         case 'decline'
%             if RLS(i) < RLS(i-1)
%                 slope = 'incline';
%             else
%                 slope = 'decline';
%             end
%     end
%     
% %     if endPtIndex == 3628
% %         scatter(rls_time(i), RLS(i));
% %     end
%     i = i-1;
% end
end

