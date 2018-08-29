% Train an SVM to classify perceptual state during rivalry

loadPaths % add dependency folders to path

clearvars

tic

% load parameters (includes SVM features)
analysisParams

% create temporary params file to use for this run
paramsFiles = dir('analysisParams*.m');
if length(paramsFiles) > 1
    lastRun = str2double(paramsFiles(end).name(15));
else
    lastRun = 0;
end
thisRun = num2str(lastRun + 1);
copyfile('analysisParams.m', ['analysisParams' thisRun '.m']);
paramsFile = ['analysisParams' thisRun '.m'];

%% Specify training and testing sets

trainSet = 'rivalry';
% 'sim' 'rivalry' 'rivalTransitions' 'simTransitions'
testSet = 'rivalTransitions';
% 'identical' 'subset' 'rivalry' 'rivalTransitions' 'simTransitions'

%% Number of states

numStates = 3;
% 2 for mixed vs. dom
% 3 for mixed vs. high dom vs. low dom

%% Specify features

% specify features in the last section of analysisParams.m

%% EDIT THIS SECTION TO CHANGE PARAMETERS

% latency
rivLat = 0; % (seconds)

% middle of DOMINANT TRAINING segment (number of seconds before end of epoch)
domSegMid_sim = 0.25; %s
domSegMid_riv = 0.4; % this is a PROP %0.5; %s

% middle of MIXED TRAINING segment (proportion of epoch length)
mixSegProp = 0.5; 

% testing timepoints within transition timecourse
binSize =  0.5; % seconds
binInt =  0.2; % interval between start of each bin (seconds)
thPoint = 1; % test on every _th time point

%% Segment parameters

segLength = 0.5; %0.5; % s
segPoints = round(segLength*sampRate); % segment length in data points (each segment is a "trial" to be classified)


%% TESTING timecourse parameters

binPts = round(binSize*sampRate); % number of points in a testing bin
intPts = round(binInt*sampRate); % number of points between start of each bin

totalTPoints = sampRate*svmTransHalf*2 + 1;
firstPt = sampRate*(transHalf-svmTransHalf) + 1;
lastPt = sampRate*(transHalf-svmTransHalf) + 1 + totalTPoints;

% set start and end points of bins
teStartPts = firstPt:intPts*thPoint:lastPt - binPts + 1;
%teStartPts = [teStartPts(15)]; % TEMP for testing
%teStartPts = 1536:1539; %[2300 2813]; %[1790 2433]; % TEMP for testing
teEndPts = teStartPts + binPts - 1;

% testing time axis
binHalf = binSize/2;
teTime = -svmTransHalf+binHalf:binInt*thPoint:svmTransHalf-binHalf;

if isempty(strfind(testSet, 'Transitions'))
    teStartPts = NaN;
    teEndPts = NaN;
end

%% Prepare run ID and directory to save data

saveDir = 'ML/groupAccMats/';

% Get the most recent run ID
files = dir([saveDir 'run*']);
if ~isempty(files)
    id = files(end).name(4:6);
else
    id = '0';
end

% Create new ID
newID = str2double(id) + 1;
runID = sprintf('%03u',newID);

%% Iterate over testing transition timepoints

accuracyMat = [];
domAccMat = [];
mixAccMat = [];
postProbMat = [];

teMsglength = 0;

for iBin = 1:length(teStartPts)
    
    teStartPt = teStartPts(iBin);
    teEndPt = teEndPts(iBin);
    
    fprintf(repmat('\b',1,teMsglength));
    teMsg = sprintf('Running testing timepoint %d of %d\n', iBin, length(teStartPts));
    fprintf(teMsg);
    teMsglength = numel(teMsg);
    
    [accuracies, postProbs, confusionMats, mixedAccs, domAccs] = svmClassify_withMixed(trainSet, testSet, rivLat, segPoints, domSegMid_sim, domSegMid_riv, mixSegProp, teStartPt, teEndPt, numStates, paramsFile);
    % trainsets: 'sim' 'rivalry' 'rivalTransitions' 'simTransitions'
    % testsets: 'identical' 'subset' 'rivalSegs' 'rivalTransitions' 'simTransitions'
    
    % store accuracies
    accuracyMat = [accuracyMat accuracies];
    domAccMat = [domAccMat domAccs];
    mixAccMat = [mixAccMat mixedAccs];
    
    % store posterior probabilities
    postProbMat = [postProbMat postProbs];
    
end

% save accuracies & posterior probabilities for each participant
save([saveDir group '_train-' trainSet '_test-' testSet '_' num2str(length(accuracies)) 'pars_3stateMixed_numPointsClassThreshold'], 'accuracyMat', 'postProbMat', 'mixAccMat', 'domAccMat', 'teTime');

delete(paramsFile);

toc