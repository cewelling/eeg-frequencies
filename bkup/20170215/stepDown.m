function [ platCount, found ] = stepDown( RLS, i, platCount, step, found)
%stepDown.m Find the closest trough
%   Detailed explanation goes here

%load parameters
analysisParams

% reached a local maximum
if RLS(i) < RLS(i+step) 
    found = 1;
    return;
    
% on a plateau
elseif RLS(i) == RLS(i+step)
    platCount = platCount + 1;
    
% still going down
else 
    platCount = 0;
end
end

