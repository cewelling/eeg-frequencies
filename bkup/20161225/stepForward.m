function [ peakIndex ] = stepForward( RLS, epochMax, i, j, forSlope, backSlope, platCount, step )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%load parameters
analysisParams

if i < length(RLS) - (timeWindow / 1000 * 512) % haven't yet reached end of the time series
    switch forSlope
        case 'incline'
            if RLS(i) >= epochMax
                if RLS(i) > RLS(i+step)
                    peakIndex = i - round(platCount/2); % choose middle of plateau
                    return;
                elseif RLS(i) == RLS(i+step)
                    platCount = platCount + 1;
                    forSlope = 'incline';
                end
            elseif RLS(i) > RLS(i+step)
                platCount = 0;
                forSlope = 'decline';
            else
                platCount = 0;
                forSlope = 'incline';
            end
            
        case 'decline'
            if RLS(i) < RLS(i+step)
                forSlope = 'incline';
            else
                forSlope = 'decline';
            end
    end
else
    % Reached end of time series; only pursue backward direction
    if j > 1
%         if j == 12215
%             disp('hi')
%         end
%         disp(num2str(j));
        peakIndex = stepBack(RLS, epochMax, i, j-step, forSlope, backSlope, platCount, -step);
        return
    else
        peakIndex = NaN;
        return
    end
    %     if isnan(peakIndex)
    %         scatter(i, RLS(i));
    %     end
end

peakIndex = stepBack(RLS, epochMax, i, j-step, forSlope, backSlope, platCount, -step);

end

