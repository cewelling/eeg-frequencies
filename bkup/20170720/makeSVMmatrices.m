function [trainingData, testData, trainingLabels, testLabels, labelVect] = makeSVMmatrices(parName, segPoints, trainSet, testSet, domSegMid_sim, domSegMid_riv, mixSegProp, teStartPt, teEndPt, rivLat, numStates, paramsFile)
% makeSVMmatrices puts EEG data into the correct format for classification
% with MatLab's SVM functions
% dependencies: getSegs.m, getSegPresses.m
% ACCOUNTS for cutting 1 second off RLS data!!!

%clearvars

%% set-up

% load parameters
run(paramsFile);

% transition bin must have at least this many zeroes to be marked as mixed
%zeroThresh = segPoints*0.3; % 30%
classThresh = segPoints*0.7; % transition bin must be at least 70% a particular state to be called that state

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

rivLat_str = num2str(round(rivLat*1000)); % for use in filenames

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
        domSegMid = domSegMid_sim;
    elseif strcmp(trainSet, 'rivalry')
        trainType = 1;
        domSegMid = domSegMid_riv;
    end
    
    trainingData = [];
    trainingLabels = [];
    for channel = channels
        for freqType = freqTypes
            
            [lowSegList, highSegList, labelVect] = getSegs(parName, date, trainType, domSegMid, mixSegProp, segPoints, channel, freqType{1}, numStates, paramsFile);
            
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
            
            if strcmp(trainingSet, 'rivalTransitions')
                if ~exist(['ML/tMatrices/' parName '_' freqType{1} '_' num2str(channel) '_' runTypes{runType} '_' rivLat_str 'msLat.mat'], 'file')
                    [lowMat, highMat, labelMat] = getTmatrix(parName, date, runType, channel, freqType, rivLat, paramsFile);
                else
                    load(['ML/tMatrices/' parName '_' freqType{1} '_' num2str(channel) '_' runTypes{runType} '_' rivLat_str 'msLat.mat'])
                end
            elseif strcmp(trainingSet, 'simTransitions')
                if ~exist(['ML/tMatrices/' parName '_' freqType{1} '_' num2str(channel) '_' runTypes{runType} '.mat'], 'file')
                    [lowMat, highMat, labelMat] = getTmatrix(parName, date, runType, channel, freqType, rivLat, paramsFile);
                else
                    load(['ML/tMatrices/' parName '_' freqType{1} '_' num2str(channel) '_' runTypes{runType} '.mat'])
                end
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
    if numStates == 3
        trainingLabels(logical(labelVect < 0)) = {'lo'};
        trainingLabels(logical(labelVect > 0)) = {'hi'};
    elseif numStates == 2
        trainingLabels(logical(labelVect ~= 0)) = {'dom'};
    end
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
if strcmp(testSet, 'rivalry')
    
    testData = [];
    testLabels = [];
    for channel = channels
        for freqType = freqTypes
            
            [lowSegList, highSegList, labelVect] = getSegs(parName, date, 1, domSegMid_riv, mixSegProp, segPoints, channel, freqType{1}, numStates, paramsFile);
            
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
            
            if strcmp(testSet, 'rivalTransitions')
                if ~exist(['ML/tMatrices/' parName '_' freqType{1} '_' num2str(channel) '_' runTypes{runType} '_' rivLat_str 'msLat.mat'], 'file')
                    [lowMat, highMat, labelMat] = getTmatrix(parName, date, runType, channel, freqType, rivLat, paramsFile);
                else
                    load(['ML/tMatrices/' parName '_' freqType{1} '_' num2str(channel) '_' runTypes{runType} '_' rivLat_str 'msLat.mat'])
                end
            elseif strcmp(testSet, 'simTransitions')
                if ~exist(['ML/tMatrices/' parName '_' freqType{1} '_' num2str(channel) '_' runTypes{runType} '.mat'], 'file')
                    [lowMat, highMat, labelMat] = getTmatrix(parName, date, runType, channel, freqType, rivLat, paramsFile);
                else
                    load(['ML/tMatrices/' parName '_' freqType{1} '_' num2str(channel) '_' runTypes{runType} '.mat'])
                end
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
            oldTestData = testData;
            if d1dt
                testData = [oldTestData lowD1Vect highD1Vect];
            end
            if d2dt
                testData = [testData lowD2Vect highD2Vect];
            end
        end
    end
    
    % test labels
    testLabelMat = labelMat(:,teStartPt:teEndPt);
    labelVect = nanmean(testLabelMat,2);
    
    % determine how many of each state each transition has (in this time bin)
    numLow = nan(size(testLabelMat, 1), 1);
    numMixed = nan(size(testLabelMat, 1), 1);
    numHigh = nan(size(testLabelMat, 1), 1);
    for i = 1:size(testLabelMat, 1)
        numLow(i) = length(find(testLabelMat(i,:) == -1));
        numMixed(i) = length(find(testLabelMat(i,:) == 0));
        numHigh(i) = length(find(testLabelMat(i,:) == 1));
    end
    
    % How many 
%     [zeroRows, zeroCols] = find(~isnan(testLabelMat) & abs(testLabelMat) - 1 ~= 0);
%     hasMixedVect = nan(length(labelVect),1);
%     hasMixedVect(hasZeroes) = 1; 
    
    % remove transition pieces with no reporting at all
    labelVect_all = labelVect;
    labelVect = labelVect(~isnan(labelVect_all));
    testData = testData(~isnan(labelVect_all),:);
    numLow = numLow(~isnan(labelVect_all));
    numMixed = numMixed(~isnan(labelVect_all));
    numHigh = numHigh(~isnan(labelVect_all));
    %hasMixedVect = hasMixedVect(~isnan(labelVect_all));
    
    testLabels = nan(size(labelVect, 1),1);
    testLabels = num2cell(testLabels);
    
    if numStates == 3
        testLabels(numLow >= classThresh) = {'lo'};
        testLabels(numMixed >= classThresh) = {'mi'};
        testLabels(numHigh >= classThresh) = {'hi'}; 
        
%         testLabels(logical(labelVect < -0.5)) = {'lo'}; %% Was < 0, > 0
%         testLabels(logical(labelVect > 0.5)) = {'hi'};
%         
%         testLabels(logical(labelVect >= -0.5 & labelVect <= 0.5 & numMixed >= zeroThresh)) = {'mi'}; %% was == 0
%         testLabels(logical(labelVect >= -0.5 & labelVect < -0.3 & numMixed < zeroThresh)) = {'lo'};
%         testLabels(logical(labelVect > 0.3 & labelVect <= 0.5 & numMixed < zeroThresh)) = {'hi'};
        
        % remove transition bins with undetermined labels (e.g. half one dom, half
        % other dom)
        %undeterm = logical(labelVect == 0 & hasMixedVect ~= 1);
        %undeterm = logical(labelVect >= -0.3 & labelVect <= 0.3 & numMixed < zeroThresh);
        undeterm = logical(numLow < classThresh & numMixed < classThresh & numHigh < classThresh);
        testLabels = testLabels(~undeterm);
        testData = testData(~undeterm,:);
    elseif numStates == 2
        testLabels(logical(labelVect ~= 0)) = {'dom'};
    end
    %     testLabels(logical(labelVect >= -0.5 & labelVect <= 0.5 & hasMixedVect == 1)) = {'mi'}; %% was == 0
    %     testLabels(logical(labelVect >= -0.5 & labelVect < 0 & hasMixedVect ~= 1)) = {'lo'};
    %     testLabels(logical(labelVect > 0 & labelVect <= 0.5 & hasMixedVect ~= 1)) = {'hi'};
    
    % limit testing data to instances where button press indicated consistent dominant state
    %testLabels = testLabels((labelVect == 1 | labelVect == -1),:);
    %testData = testData((labelVect == 1 | labelVect == -1),:);
    
     % remove transition bins that are all NaNs (because they intersected with
    % the beginning or end of a run)
    testData_all = testData;
    testData = testData(~isnan(nanmean(testData_all, 2)),:);
    testLabels = testLabels(~isnan(nanmean(testData_all, 2)),:);
end
end