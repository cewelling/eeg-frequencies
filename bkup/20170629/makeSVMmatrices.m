function [trainingData, testData, trainingLabels, testLabels] = makeSVMmatrices(parName, segPoints, trainSet, testSet, segMiddle, trStartPt, trEndPt, teStartPt, teEndPt, rivLat)
% makeSVMmatrices puts EEG data into the correct format for classification
% with MatLab's SVM functions
% dependencies: getSegs.m, getSegPresses.m
% ACCOUNTS for cutting 1 second off RLS data!!!

%clearvars

%% set-up

% load parameters
analysisParams

% Handle inputs for timecourse options
if strfind(trainSet, 'Transitions')
    if isnan(trStartPt)
        error('Provide a startpt for the training bin');
    elseif isnan(trEndPt)
        error('Provide an endpt for the training bin');
    end
elseif strfind(testSet, 'Transitions')
    if isnan(teStartPt)
        error('Provide a startpt for the testing bin');
    elseif isnan(teEndPt)
        error('Provide an endpt for the testing bin');
    end
end

rivLat_str = num2str(rivLat); % for use in filenames

% Get EEG file info
EEGfiles = dir([eegDir parName '*']);
% ...if the participant exists...
if isempty(EEGfiles)
    return;
end
date = strtok(EEGfiles(1).date);
date = datestr(date, dateFormat);

%% Create training set: segments
if strcmp(trainSet, 'sim') || strcmp(trainSet, 'rivalry')
    
    if strcmp(trainSet, 'sim')
        trainType = 2;
    elseif strcmp(trainSet, 'rivalry')
        trainType = 1;
    end
    
    trainingData = [];
    trainingLabels = [];
    for channel = channels
        for freqType = freqTypes
            
            [lowSegList, highSegList, labelVect] = getSegs(parName, date, trainType, segMiddle, segPoints, channel, freqType{1});
            
            % Bin the data according to the given bin interval
            %                 leftOver = mod(segPoints, binPoints);
            %                 for currPoint = ceil(leftOver/2)+binPoints:binPoints:segPoints - floor(leftOver/2)
            %                     lowBinVal = nanmean(lowSeg(currPoint - binPoints:currPoint - 1));
            %                     highBinVal = nanmean(highSeg(currPoint - binPoints:currPoint - 1));
            %                     lowFreqVect = [lowFreqVect; lowBinVal];
            %                     highFreqVect = [highFreqVect; highBinVal];
            %                 end
            
            % average over segments, calculate derivatives
            if ~isempty(labelVect)
                %plotSegs;
                lowFreqVect = nanmean(lowSegList, 2);
                lowD1Vect = nanmean(diff(lowSegList, 1, 2), 2);
                lowD2Vect = nanmean(diff(lowSegList, 2, 2), 2);
                highFreqVect = nanmean(highSegList, 2);
                highD1Vect = nanmean(diff(highSegList, 1, 2), 2);
                highD2Vect = nanmean(diff(highSegList, 2, 2), 2);
                
                if amp
                    trainingData = [trainingData lowFreqVect highFreqVect];
                end
                if d1dt
                    trainingData = [trainingData lowD1Vect highD1Vect];
                end
                if d2dt
                    trainingData = [trainingData lowD2Vect highD2Vect];
                end
                
                trainingLabels = labelVect; % labelVect should be the same for all frequencies / channels
            end
        end
    end
end

%% Create training set: transitions
if strcmp(trainSet, 'rivalTransitions') || strcmp(trainSet, 'simTransitions')
    
    if strcmp(trainSet, 'rivalTransitions')
        runType = 1; % use rivalry transitions
    else
        runType = 2; % use simulation transitions
    end
    
    trainingData = [];
    trainingLabels = [];
    
    for channel = channels
        for freqType = freqTypes
            
            if ~exist(['ML/tMatrices/' parName '_' freqType{1} '_' num2str(channel) '_' runTypes{runType} '_' rivLat_str(3:end) 'msLat.mat'], 'file')
                [lowMat, highMat, labelMat] = getTmatrix(parName, date, runType, channel, freqType, rivLat);
            else
                load(['ML/tMatrices/' parName '_' freqType{1} '_' num2str(channel) '_' runTypes{runType} '_' rivLat_str(3:end) 'msLat.mat'])  
            end
            
            % training data
            lowFreqVect = nanmean(lowMat(:,trStartPt:trEndPt),2);
            lowD1Vect = nanmean(diff(lowMat(:,trStartPt:trEndPt), 1, 2),2);
            lowD2Vect = nanmean(diff(lowMat(:,trStartPt:trEndPt), 2, 2),2);
            highFreqVect = nanmean(highMat(:,trStartPt:trEndPt),2);
            highD1Vect = nanmean(diff(highMat(:,trStartPt:trEndPt), 1, 2),2);
            highD2Vect = nanmean(diff(highMat(:,trStartPt:trEndPt), 2, 2),2);
            
            if amp
                trainingData = [trainingData lowFreqVect highFreqVect];
            end
            if d1dt
                trainingData = [trainingData lowD1Vect highD1Vect];
            end
            if d2dt
                trainingData = [trainingData lowD2Vect highD2Vect];
            end
        end
    end
    
    % training labels
    labelVect = nanmean(labelMat(:,trStartPt:trEndPt),2);
    trainingLabels = cell(size(labelVect,1),1);
    trainingLabels(logical(labelVect < 0)) = {'lo'};
    trainingLabels(logical(labelVect > 0)) = {'hi'};
    trainingLabels(logical(labelVect == 0)) = {'mi'};
    
    % limit training data to instances where button press indicated consistent dominant state
    %trainingLabels = trainingLabels((labelVect == 1 | labelVect == -1),:);
    %trainingData = trainingData((labelVect == 1 | labelVect == -1),:);
    
    % remove transitions that are all NaNs (because they intersected with
    % the beginning or end of a run)
    trainingData_all = trainingData;
    trainingData = trainingData(~isnan(nanmean(trainingData_all, 2)),:);
    trainingLabels = trainingLabels(~isnan(nanmean(trainingData_all, 2)),:);
end

%% Create test set (identical to training set)
if strcmp(testSet, 'identical')
    testData = trainingData;
    testLabels = trainingLabels;
end

%% Create test set (subset of sim or rivalry trials)
if strcmp(testSet, 'subset')
    splits = 5; % 4/5 for training, 1/5 for testing
    
    % set indices for training and test data
    
    % random method
    %     totalInds = 1:length(trainingLabels);
    %     testInds = randperm(length(trainingLabels), round(length(trainingLabels)/splits));
    %     trainingInds = setdiff(totalInds, testInds);
    
    % set method
    totalInds = 1:length(trainingLabels);
    scrambled = randperm(length(trainingLabels));
    trainingInds = totalInds(mod(scrambled, splits) ~= 0);
    testInds = totalInds(mod(totalInds, splits) == 0);
    
    % split data
    testData = trainingData(testInds,:);
    testLabels = trainingLabels(testInds,:);
    trainingData = trainingData(trainingInds,:);
    trainingLabels = trainingLabels(trainingInds,:);
end

%% Create test set (rivalry segments)
if strcmp(testSet, 'rivalSegs')
    segMiddle = 0.6; % middle of the segment in seconds after epoch start
    
    testData = [];
    testLabels = [];
    for channel = channels
        for freqType = freqTypes
            
            [lowSegList, highSegList, labelVect] = getSegs(parName, date, 1, segMiddle, segPoints, channel, freqType{1});
            
            % Bin the data according to the given bin interval
            %                 leftOver = mod(segPoints, binPoints);
            %                 for currPoint = ceil(leftOver/2)+binPoints:binPoints:segPoints - floor(leftOver/2)
            %                     lowBinVal = nanmean(lowSeg(currPoint - binPoints:currPoint - 1));
            %                     highBinVal = nanmean(highSeg(currPoint - binPoints:currPoint - 1));
            %                     lowFreqVect = [lowFreqVect; lowBinVal];
            %                     highFreqVect = [highFreqVect; highBinVal];
            %                 end
            
            if ~isempty(labelVect)
                %plotSegs;
                % average over segments, calculate derivatives
                lowFreqVect = nanmean(lowSegList, 2);
                lowD1Vect = nanmean(diff(lowSegList, 1, 2), 2);
                lowD2Vect = nanmean(diff(lowSegList, 2, 2), 2);
                highFreqVect = nanmean(highSegList, 2);
                highD1Vect = nanmean(diff(highSegList, 1, 2), 2);
                highD2Vect = nanmean(diff(highSegList, 2, 2), 2);
                
                if amp
                    testData = [testData lowFreqVect highFreqVect];
                end
                if d1dt
                    testData = [testData lowD1Vect highD1Vect];
                end
                if d2dt
                    testData = [testData lowD2Vect highD2Vect];
                end
                
                testLabels = labelVect;
            end
        end
    end
end

%% Create test set (one timepoint in transitions)
if strcmp(testSet, 'rivalTransitions') || strcmp(testSet, 'simTransitions')
    
    if strcmp(testSet, 'rivalTransitions')
        runType = 1; % use rivalry transitions
    else
        runType = 2; % use simulation transitions
    end
    
    testData = [];
    testLabels = [];
    
    for channel = channels
        for freqType = freqTypes
            
            if ~exist(['ML/tMatrices/' parName '_' freqType{1} '_' num2str(channel) '_' '.mat'], 'file')
                [lowMat, highMat, labelMat] = getTmatrix(parName, date, runType, channel, freqType, rivLat);
            else
                load(['ML/tMatrices/' parName '_' freqType{1} '_' num2str(channel) '.mat'])
            end
            
            % test data
            lowFreqVect = nanmean(lowMat(:,teStartPt:teEndPt),2);
            lowD1Vect = nanmean(diff(lowMat(:,teStartPt:teEndPt), 1, 2),2);
            lowD2Vect = nanmean(diff(lowMat(:,teStartPt:teEndPt), 2, 2),2);
            highFreqVect = nanmean(highMat(:,teStartPt:teEndPt),2);
            highD1Vect = nanmean(diff(highMat(:,teStartPt:teEndPt), 1, 2),2);
            highD2Vect = nanmean(diff(highMat(:,teStartPt:teEndPt), 2, 2),2);
            
            if amp
                testData = [testData lowFreqVect highFreqVect];
            end
            if d1dt
                testData = [testData lowD1Vect highD1Vect];
            end
            if d2dt
                testData = [testData lowD2Vect highD2Vect];
            end
        end
    end
    
    % test labels
    labelVect = nanmean(labelMat(:,teStartPt:teEndPt),2);
    testLabels = cell(size(labelVect,1),1);
    testLabels(logical(labelVect < 0)) = {'lo'};
    testLabels(logical(labelVect > 0)) = {'hi'};
    testLabels(logical(labelVect == 0)) = {'mi'};
    
    % limit testing data to instances where button press indicated consistent dominant state
    %testLabels = testLabels((labelVect == 1 | labelVect == -1),:);
    %testData = testData((labelVect == 1 | labelVect == -1),:);
    
     % remove transitions that are all NaNs (because they intersected with
    % the beginning or end of a run)
    testData_all = testData;
    testData = testData(~isnan(nanmean(testData_all, 2)),:);
    testLabels = testLabels(~isnan(nanmean(testData_all, 2)),:);
end
end