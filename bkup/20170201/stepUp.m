function [ platCount, found ] = stepUp( RLS, i, platCount, step, found)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%load parameters
analysisParams

% reached a local maximum
if RLS(i) > RLS(i+step) 
    found = 1;
    return;
    
% on a plateau
elseif RLS(i) == RLS(i+step)
    platCount = platCount + 1;
    
% still going up
else 
    platCount = 0;
end
end

