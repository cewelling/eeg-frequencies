function [ newStruct ] = structSubset( structure, ind)
%UNTITLED5 Summary of this function goes here
%   Takes a 1x1 transition list structure, where each field is a matrix 
% with the same, number of rows and returns the subset of that structure
% (the rows from ind to the end). 
% 
% Structures must have the following fields:
% f1
% f2
% tPoints
% tTrials

if size(structure) ~= [1 1]
    error('Input structure must be 1x1')
end

newStruct.f1 = structure.f1(ind:end,:);
newStruct.f2 = structure.f2(ind:end,:);
newStruct.tPoints = structure.tPoints(ind:end,:);
newStruct.tTrials = structure.tTrials(ind:end,:);

