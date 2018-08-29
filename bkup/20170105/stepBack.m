function [ backSlope, platCount, found ] = stepBack(RLS, epochMax, j, backSlope, platCount, found)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%load parameters
analysisParams

switch backSlope
    case 'incline'
        % reached a point higher than absolute max of epoch
        if RLS(j) >= epochMax
            if RLS(j) > RLS(j-1)
                found = 1; % choose middle of plateau
                return;
            elseif RLS(j) == RLS(j-1)
                platCount = platCount + 1;
                backSlope = 'incline';
            end
        % reached a local maximum, but not the highest of this epoch
        elseif RLS(j) > RLS(j-1)
            platCount = 0;
            backSlope = 'decline';
        % still going up
        else
            platCount = 0;
            backSlope = 'incline';
        end
        
    case 'decline'
        if RLS(j) < RLS(j-1)
            backSlope = 'incline';
        else
            backSlope = 'decline';
        end
end
end

