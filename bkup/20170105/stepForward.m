function [ forSlope, platCount, found ] = stepForward( RLS, epochMax, i, forSlope, platCount, found)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%load parameters
analysisParams

switch forSlope
    case 'incline'
        % reached a point higher than absolute max of epoch
        if RLS(i) >= epochMax
            if RLS(i) > RLS(i+1)
                found = 1;
                return;
            elseif RLS(i) == RLS(i+1)
                platCount = platCount + 1;
                forSlope = 'incline';
            end
        % reached a local maximum, but not the highest of this epoch
        elseif RLS(i) > RLS(i+1)
            platCount = 0;
            forSlope = 'decline';
        % still going up
        else
            platCount = 0;
            forSlope = 'incline';
        end
        
    case 'decline'
        if RLS(i) < RLS(i+1)
            forSlope = 'incline';
        else
            forSlope = 'decline';
        end
end
end

