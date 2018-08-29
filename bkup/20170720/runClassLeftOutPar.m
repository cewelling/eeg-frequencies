% Train an SVM to classify dominant state of a left out participant using
% button presses of all other participants 
% Should this work across groups? Yes I think...
% Dependencies: svmClassLeftOutPar

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
testSet = 'identical'; % identical: same as trainSet, but in the left out participant
% 'identical' 'subset' 'rivalry' 'rivalTransitions' 'simTransitions' 
% 'rivalry': rivalry segments

%% Run with mixed states?
withMixed = 1;
numStates = 3;
% 2 for mixed vs. dom
% 3 for mixed vs. high dom vs. low dom (so if classifying high dom vs. low dom, pick this option)

%% Specify features

% specify features in the last section of analysisParams.m

%% EDIT THIS SECTION TO CHANGE PARAMETERS

% latency
rivLat = 0; % (seconds)

% middle of DOMINANT TRAINING segment (number of seconds before end of epoch)
domSegMid_sim = 0.25; %s
domSegMid_riv = 0.4; % this is actually segProp %0.5; %s

% middle of MIXED TRAINING segment (proportion of epoch length)
mixSegProp = 0.5; 

% testing timepoints within transition timecourse
binSize = 0.5; % seconds
binInt = 0.2; % interval between start of each bin (seconds)
thPoint = 1; % test on every _th time point

%% Segment parameters

segLength = 0.5; %0.5; % s
segPoints = round(segLength*sampRate); % segment length in data points (each segment is a "trial" to be classified)


%% TESTING timecourse parameters

if ~isempty(strfind(testSet, 'Transitions'))
    binPts = round(binSize*sampRate); % number of points in a testing bin
    intPts = round(binInt*sampRate); % number of points between start of each bin
    
    totalTPoints = sampRate*svmTransHalf*2 + 1;
    firstPt = sampRate*(transHalf-svmTransHalf) + 1;
    lastPt = sampRate*(transHalf-svmTransHalf) + 1 + totalTPoints;
    
    % set start and end points of bins
    teStartPts = firstPt:intPts*thPoint:lastPt - binPts + 1;
    %teStartPts = 2298; %[2300 2813]; %[1790 2433]; % TEMP for testing
    teEndPts = teStartPts + binPts - 1;
    
    % testing time axis
    binHalf = binSize/2;
    teTime = -svmTransHalf+binHalf:binInt*thPoint:svmTransHalf-binHalf;
    
else
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
    
    [accuracies, postProbs, confusionMats, mixedAccs, domAccs] = svmClassLeftOutPar(trainSet, testSet, rivLat, segPoints, domSegMid_sim, domSegMid_riv, mixSegProp, teStartPt, teEndPt, paramsFile, withMixed, numStates);
    
    % trainsets: 'sim' 'rivalry' 'rivalTransitions' 'simTransitions'
    % testsets: 'identical' 'subset' 'rivalSegs' 'rivalTransitions' 'simTransitions'
    
    if ~isempty(strfind(testSet, 'Transitions'))
        
        % store accuracies
        accuracyMat = [accuracyMat accuracies];
        domAccMat = [domAccMat domAccs];
        mixAccMat = [mixAccMat mixedAccs];
        
        % store posterior probabilities
        postProbMat = [postProbMat postProbs]; 
    else
        
        % plot histogram
        disp(['group: ' group])
        disp(['mean: ' num2str(nanmean(accuracies))])
        disp(['median: ' num2str(nanmedian(accuracies))])
        disp(['min: ' num2str(nanmin(accuracies))])
        disp(['max: ' num2str(nanmax(accuracies))])
        disp(['']);
        
        figure
        histogram(accuracies, 15);
        xlabel('accuracy')
        title('Accuracies: classify a left out participant')
    end
end

% save accuracies & posterior probabilities for each participant
if ~isempty(strfind(testSet, 'Transitions'))
    if withMixed
        save([saveDir 'bothGroups_train-' trainSet '_test-' testSet '_' num2str(length(accuracies)) 'left_out_pars_withMixed'], 'accuracyMat', 'postProbMat', 'domAccMat', 'mixAccMat', 'teTime');
    else
        save([saveDir 'bothGroups_train-' trainSet '_test-' testSet '_' num2str(length(accuracies)) 'left_out_pars_PROP'], 'accuracyMat', 'postProbMat', 'domAccMat', 'mixAccMat', 'teTime');
    end
end

delete(paramsFile)

toc
