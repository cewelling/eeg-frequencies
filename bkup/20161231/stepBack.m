function [ peakIndex ] = stepBack(RLS, epochMax, i, j, forSlope, backSlope, platCount, step)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%load parameters
analysisParams

if j > 1 % haven't yet reached beginning of the time series
%     disp(num2str(j));
    switch backSlope
        case 'incline'
            if RLS(j) >= epochMax
                if RLS(j) > RLS(j+step)
                    peakIndex = j - round(platCount/2); % choose middle of plateau
                    return;
                elseif RLS(j) == RLS(j+step)
                    platCount = platCount + 1;
                    backSlope = 'incline';
                end
            elseif RLS(j) > RLS(j+step)
                platCount = 0;
                backSlope = 'decline';
            else
                platCount = 0;
                backSlope = 'incline';
            end
            
        case 'decline'
            if RLS(j) < RLS(j+step)
                backSlope = 'incline';
            else
                backSlope = 'decline';
            end
    end
    
else
    % beginning of time series reached; only pursue forward direction
    if i < length(RLS) - (timeWindow / 1000 * 512)
        peakIndex = stepForward(RLS, epochMax, i-step, j, forSlope, backSlope, platCount, -step);
        return
    else
        peakIndex = NaN;
        return
    end
    %     if isnan(peakIndex)
    %         scatter(i, RLS(i));
    %     end
end

%if i < length(RLS) - (timeWindow / 1000 * 512) .....
peakIndex = stepForward(RLS, epochMax, i-step, j, forSlope, backSlope, platCount, -step);

end

