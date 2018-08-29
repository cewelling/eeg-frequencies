%% Parameters for SVM scripts

%% General Params
acrossPars = 0; %Classify across or within participant? 1 for across, 0 for within 
classes = '3state'; % 'hiVlo' % 'mixVdom' % '3state'

trainSet = 'rivalry'; % 'sim' 'rivalry' 'rivalTransitions' 'simTransitions'
testSet = 'identical'; % 'identical' 'subset' 'rivalry' 'rivalTransitions' 'simTransitions' 
%subset means 1/5ths of the data specified in train
%rivalry means rivalrySegs

binSize = 0.5; % seconds
binInt = 0.2; % interval between start of each bin (seconds)
thPoint = 1; % test on every _th time point

svmTransHalf = 3; % how many second to the L/R of the middle of the transition do we test/plot

minTrEx = 30; % minimum number of examples in training set for inclusion

%% Features for SVM
amp = 0; % amplitude
d1dt = 0; % 1st derivative
d2dt = 1; % 2nd derivative
freqTypes = {'principles'}; %{'principles' 'harmonics' 'im'}; % 'harmonics' 'im'}; %{'principles' 'harmonics' 'im'};
channels = focusElectrode.nums; %occipitals.nums; %focusElectrode.nums;

%% Pretty locked down General Params 
% latency
rivLat = 0; % (seconds)

% middle of DOMINANT TRAINING segment (number of seconds before end of epoch)
domSegMid_sim = 0.25; %s
domSegMid_riv = 0.4; % this is actually segProp %0.5; %s

% middle of MIXED TRAINING segment (proportion of epoch length)
mixSegProp = 0.5; 

minEpochDur = 0.5; % minimum duration to be considered an epoch (s)

segLength = 0.5; %0.5; % s
segPoints = round(segLength*sampRate); % segment length in data points 


%% Switches
% feature visualization
plot2features = 1; % You must go into the script (e.g. svmClassLeftOutPar) to choose which 2 features to plot
segPlotOrNot = 'no'; % visualize segments laid over RLS traces?


%% Type of classification
if strcmp(classes, 'hiVlo')
    withMixed = 0;
    numStates = 3;
    boxConstraint = 0.01; %This was optimized using runSVMCrossVal on the LOPs?
elseif strcmp(classes, 'mixVdom')
    withMixed = 1;
    numStates = 2;
    boxConstraint = 0.01; %This was optimized using runSVMCrossVal on the LOPs?
elseif strcmp(classes, '3state')
    withMixed = 1;
    numStates = 3;
    boxConstraint = 0.01; %This was optimized using runSVMCrossVal on the LOPs?
end

%% Select participants
if ~acrossPars
    % For ****WITHIN PARTICIPANT**** analysis
    group = 0; %0 controls, 1 for asc
    if group == 0
        parNums = [105:140]; % ex/ Stratus100's parNum is 100 %Note, CER = 137
        group = 'stratus'; % 'stratus' or 'cumulus'
        % exclude participants
        load('pre-processing/stratExcluded.mat')
        parNums = setdiff(parNums,stratExcluded);
    else
        parNums = [1:28]; % ex/ Stratus100's parNum is 100
        group = 'cumulus'; % 'stratus' or 'cumulus'
        % exclude participants
        load('pre-processing/cumExcluded.mat')
        parNums = setdiff(parNums,cumExcluded);
    end
    
else
    % For ****ACROSS PARTICIPANT**** analysis
    groupParNums = {[105:138], [1:28]};
    groupCodes = {'stratus','cumulus'}; %{'stratus','cumulus'};
    analGroupIDs = {'stratus', 'cumulus'}; %{'stratus', 'cumulus'};
    
    % exclude participants
    load('pre-processing/stratExcluded.mat')
    load('pre-processing/cumExcluded.mat')
    groupParNums{1} = setdiff(groupParNums{1},stratExcluded);
    groupParNums{2} = setdiff(groupParNums{2},cumExcluded);
    group = 'bothGroups';   
end


%% TESTING timecourse parameters
if ~isempty(strfind(testSet, 'Transitions'))
    binPts = round(binSize*sampRate); % number of points in a testing bin
    intPts = round(binInt*sampRate); % number of points between start of each bin
    
    totalTPoints = sampRate*svmTransHalf*2 + 1;
    firstPt = sampRate*(transHalf-svmTransHalf) + 1;
    lastPt = sampRate*(transHalf-svmTransHalf) + 1 + totalTPoints;
    
    % set start and end points of bins
    teStartPts = firstPt:intPts*thPoint:lastPt - binPts + 1;
    %teStartPts = [2300 2813]; %[1790 2433]; % TEMP for testing
    teEndPts = teStartPts + binPts - 1;
    
    % testing time axis
    binHalf = binSize/2;
    teTime = -svmTransHalf+binHalf:binInt*thPoint:svmTransHalf-binHalf;
    
else
    teStartPts = NaN;
    teEndPts = NaN;
end

% use sim schedule (rather than button presses) for sim trials? 
% (for both transitions and segments)
useSimSchedule = 1;  