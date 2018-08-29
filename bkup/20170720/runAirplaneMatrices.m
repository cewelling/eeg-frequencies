% optimize training segment and latency for classification of dominant
% rivalry state with an SVM

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
% 'rivalry': rivalry segments

%% Run with mixed states?
withMixed = 0;

numStates = 3;
% 2 for mixed vs. dom
% 3 for mixed vs. high dom vs. low dom

%% Specify features

% specify features in the last section of analysisParams.m

%% Iterating variables (EDIT THIS SECTION TO CHANGE ITERATION PARAMS)

% latencies
runOneLat = 1; % one latency at a time (so we can submit a different job to openMind for each latency)
rivLats = 0; %0:0.1:0.6; % (seconds) %0.5666;
if runOneLat == 1 && size(rivLats,1) > 1
    error('Only input one rivalry latency when runOneLat is set to 1');
end

% Vary DOMINANT or MIXED training segment position
segType = 'dom'; % 'dom' % 'mix'

% middle of DOMINANT TRAINING segment (number of seconds before end of epoch)

domSegMids_sim = 0.25; %0.25:0.05:0.75; %(s) 
domSegMids_riv = 0.1:0.1:0.9; %0.5;

% middle of MIXED TRAINING segment: proportion of epoch
mixSegProps = 0.5; %0.1:0.1:0.9; 

% testing timepoints within transition timecourse
binSize = 0.5; % seconds
binInt = 0.2; % interval between start of each bin (seconds)
thPoint = 1; % test on every _th time point

%% Segment parameters

% Segment length
segLength = 0.5; %0.5; % s
segPoints = round(segLength*sampRate); % segment length in data points (each segment is a "trial" to be classified)

% Check that we've varied the appropriate type of segment, set segment
% position axis
if strcmp(segType, 'dom')
    if strcmp(trainSet, 'sim')
        if length(domSegMids_riv) > 1 || length(mixSegProps) > 1
            error(['Only vary the position of dominant training segments: ' trainSet]);
        end
        
        % set segment postion axis
        segTime = domSegMids_sim;
        
        % specify all segment positions for each iteration
        domSegMids_riv = repmat(domSegMids_riv, 1, length(segTime));
        mixSegProps = repmat(mixSegProps, 1, length(segTime));
    elseif strcmp(trainSet, 'rivalry')
        if length(domSegMids_sim) > 1 || length(mixSegProps) > 1
            error(['Only vary the position of dominant training segments: ' trainSet]);
        end
        segTime = domSegMids_riv;
        domSegMids_sim = repmat(domSegMids_sim, 1, length(segTime));
        mixSegProps = repmat(mixSegProps, 1, length(segTime));
    end
elseif strcmp(segType, 'mix')
    if length(domSegMids_riv) > 1 || length(domSegMids_sim) > 1
        error('Only vary the position of mixed training segments');
    end
    segTime = mixSegProps;
    domSegMids_sim = repmat(domSegMids_sim, 1, length(segTime));
    domSegMids_riv = repmat(domSegMids_riv, 1, length(segTime));
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
    %teStartPts = 2298; %[2300 2813]; %[1790 2433]; % TEMP for testing
    teEndPts = teStartPts + binPts - 1;
    
    % testing time axis
    binHalf = binSize/2;
    teTime = -svmTransHalf+binHalf:binInt*thPoint:svmTransHalf-binHalf;
else
    teStartPts = NaN;
    teEndPts = NaN;
    teTime = NaN;
end

%% Prepare run ID and directory to save data

% Get the most recent airplaneMats run ID
folders = dir('ML/accMats/run*');
if ~isempty(folders)
    id = folders(end).name(4:6);
else
    id = '0';
end

% Create new ID
newID = str2double(id) + 1;
runID = sprintf('%03u',newID);

% Create folder to store accuracies for individual latencies and training segments 
mkdir(['ML/accMats/run' runID]);
accMatDir = ['ML/accMats/run' runID '/'];

%% Iterate

airplaneMat = nan(length(rivLats), length(segTime), length(teStartPts));
airplaneMatMixed = nan(length(rivLats), length(segTime), length(teStartPts));

latMsgLength = 0;
trMsglength = 0;
teMsglength = 0;

% over latencies
for iLat = 1:length(rivLats)
    
    rivLat = rivLats(iLat);
    
    fprintf(repmat('\b',1,latMsgLength + trMsglength + teMsglength));
    msg = sprintf('Running latency %d of %d\n', iLat, length(rivLats));
    fprintf(msg);
    latMsglength = numel(msg);
    
    trMsglength = 0;
    teMsglength = 0;
    
    % over training segments
    for iSeg = 1:length(segTime) 
        
        domSegMid_sim = domSegMids_sim(iSeg);
        domSegMid_riv = domSegMids_riv(iSeg);
        mixSegProp = mixSegProps(iSeg);
        
        fprintf(repmat('\b',1,trMsglength + teMsglength));
        msg = sprintf('Running training segment %d of %d\n', iSeg, length(segTime));
        fprintf(msg);
        trMsglength = numel(msg);
        
        % over testing transition timepoints
        accuracyMat = []; 
        accuracyMatMixed = [];
        
        teMsglength = 0;
        
        for iBin = 1:length(teStartPts)
            
            teStartPt = teStartPts(iBin);
            teEndPt = teEndPts(iBin);
            
            fprintf(repmat('\b',1,teMsglength));
            teMsg = sprintf('Running testing timepoint %d of %d\n', iBin, length(teStartPts));
            fprintf(teMsg);
            teMsglength = numel(teMsg);
            
            if withMixed
                [accuracies, ~, confusionMats] = svmClassify_withMixed(trainSet, testSet, rivLat, segPoints, domSegMid_sim, domSegMid_riv, mixSegProp, teStartPt, teEndPt, numStates, paramsFile);
            else
                [accuracies, ~] = svmClassify(trainSet, testSet, rivLat, segPoints, domSegMid_sim, domSegMid_riv, mixSegProp, teStartPt, teEndPt, numStates, paramsFile);
                % still has a mixSegProp parameter since getSegs.m will
                % pick mixed segments, but won't use them.
            end
            
            % trainsets: 'sim' 'rivalry' 'rivalTransitions' 'simTransitions'
            % testsets: 'identical' 'subset' 'rivalSegs' 'rivalTransitions' 'simTransitions'
            
            % store accuracies
            accuracyMat = [accuracyMat accuracies];
            airplaneMat(iLat, iSeg, iBin) = nanmean(accuracies);
            
            if withMixed
                % store *mixed* accuracies
                parMixedAccs = [];
                for i = 1:length(confusionMats)
                    thisMat = confusionMats{i};
                    thisAcc = thisMat(2,2)/sum(thisMat(2,:));
                    parMixedAccs = [parMixedAccs; thisAcc];
                end
                accuracyMatMixed = [accuracyMatMixed parMixedAccs];
                airplaneMatMixed(iLat, iSeg, iBin) = nanmean(parMixedAccs);
            end
        end
        
        % save accuracies for each participant
        save([accMatDir group '_train-' trainSet '_test-' testSet '_lat' num2str(iLat) '_seg' num2str(iSeg)], 'accuracyMat', 'teTime', 'accuracyMatMixed');
    end
end



% save airplane matrix
if runOneLat
    save(['ML/airplaneMats/run' runID '_' group '_' num2str(length(accuracies)) 'pars_train-' trainSet '_test-' testSet '_' num2str(round(1000*rivLat)) 'msLat_domSegPROP'], 'airplaneMat', 'rivLats', 'segTime','teTime', 'airplaneMatMixed')
else
    save(['ML/airplaneMats/run' runID '_' group '_' num2str(length(accuracies)) 'pars_train-' trainSet '_test-' testSet], 'airplaneMat', 'rivLats', 'segTime','teTime', 'airplaneMatMixed')
end

delete(paramsFile);

toc
