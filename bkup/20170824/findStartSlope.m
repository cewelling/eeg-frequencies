function [ slope ] = findStartSlope(RLS, index, step)
%findStartSlope: At the given point in the time-frequency spectrum, was the 
% slope on the incline or decline (in the designated step direction)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Called from: findPeak.m, itself
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load parameters
analysisParams
% reached beginning of time-frequency spectrum
if step == 1 && index == 1
    slope = 'decline';
    return;
% reached end of time-frequency spectrum
elseif step == -1 && index >= length(RLS) - (cutoffWindow / 1000 * 512) 
    %end is padded to eliminate initial spike; 512 = sampling rate
    slope = 'decline';
    return;  
elseif RLS(index) == RLS(index-step)
    slope = findStartSlope(RLS, index-step);
    return;
elseif RLS(index) > RLS(index-step)
    slope = 'incline';
    return;
elseif RLS(index) < RLS(index-step)
    slope = 'decline';
    return;
end

