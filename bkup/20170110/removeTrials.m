function [ all_data ] = removeTrials( all_data, numTrials, varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% total # of trials recorded in EEG data file
totalTrials = numTrials + nargin - 2;

% specify trials to be included
trialsIn = {};
for i = 1:nargin - 2
    if i == 1 % first trial for removal
        trialsIn = {1:varargin{i}-1};
        if i == nargin - 2 % just one trial to remove
            trialsIn = [trialsIn varargin{i} + 1:totalTrials];
        end
    elseif i == nargin - 2 % last trial for removal
        trialsIn = [trialsIn varargin{i-1} + 1:varargin{i} - 1];
        trialsIn = [trialsIn varargin{i}+1:totalTrials];
    else % trials for removal in between
        trialsIn = [trialsIn varargin{i-1} + 1:varargin{i} - 1];
    end

end

% keep only trials to be included
newTrial = {};
newTime = {};
newTrialInfo = [];
newSampleInfo = [];
for i = 1:length(trialsIn)
    newTrial = [newTrial all_data.trial(trialsIn{i})];
    newTime = [newTime all_data.time(trialsIn{i})];
    newTrialInfo = [newTrialInfo ; all_data.trialinfo(trialsIn{i})];
    newSampleInfo = [newSampleInfo; all_data.sampleinfo(trialsIn{i},:)];
end

all_data.trial = newTrial;
all_data.time = newTime;
all_data.trialinfo = newTrialInfo;
all_data.sampleinfo = newSampleInfo;

end

