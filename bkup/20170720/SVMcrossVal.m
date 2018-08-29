%function [accuracies, postProbs, confusionMats] = SVMcrossVal(trainSet, testSet, rivLat, segPoints, domSegMid_sim, domSegMid_riv, mixSegProp, teStartPt, teEndPt, paramsFile, withMixed, numStates)

% Uses classifier to...
% uses Matlab's SVM toolbox

%loadPaths % add dependency folders to path

%clearvars

%tic

% load parameters
groupAnalysisParams
analysisParams;

%% SVM parameters

trainSet = 'rivalry';
testSet = 'identical';
rivLat = 0;
segPoints = 0.5 * sampRate;
domSegMid_sim = 0.25;
domSegMid_riv = 0.5;
mixSegProp = 0.5;
teStartPt = NaN;
teEndPt = NaN;
withMixed = 1;
numStates = 3; % 2

%% C and gamma values

logCs = -5:5;
Cs = 10.^(logCs); % box constraints
logSigmas = -5:5;
sigmas = 10.^(logSigmas); % for rbf model only

%% Models

models = {'linear'}; %{'linear', 'rbf'}

%% Iterate through participants to create training and testing sets

group = 'stratus';
parNums =  [106 112 131 132]; % Excluded participants for cross-validation

groupTrainingData = [];
groupTrainingLabels = [];
groupTestData = [];
groupTestLabels = [];
trainParVect = [];
testParVect = [];

for iGroup = 1:length(groupParNums)
    group = groupCodes{iGroup};
    for parNum = groupParNums{iGroup}
        
        % participant-specific set-up
        numFormat = '%02d';
        parName = [group num2str(parNum, numFormat)];
        
        % Get EEG file info to determine whether participant exists
        EEGfiles = dir([eegDir parName '*']);
        % ...if the participant exists...
        if isempty(EEGfiles)
            continue;
        end
        
        % Have training and testing data been saved for this participant?
        trainExist = exist(['ML/parTrainSets/' parName '_' trainSet '_' num2str(numStates) 'states.mat'], 'file');
        if trainExist
            load(['ML/parTrainSets/' parName '_' trainSet '_' num2str(numStates) 'states.mat'])
        end
        if strcmp(testSet, 'identical')
            testExist = exist(['ML/parTestSets/' parName '_' trainSet '_' num2str(numStates) 'states.mat'], 'file');
            if testExist
                load(['ML/parTestSets/' parName '_' trainSet '_' num2str(numStates) 'states.mat']);
            end
        elseif strcmp(testSet, 'subset')
            error('Why are you running left out participant classification with a data subset?');
        else
            testExist = exist(['ML/parTestSets/' parName '_' testSet '_' num2str(numStates) 'states.mat'], 'file');
            if testExist
                load(['ML/parTestSets/' parName '_' testSet '_' n2str(numStates) 'states.mat']);
            end
        end
        
        if ~(trainExist && testExist)
            % Get training and testing data for this participant
            [trainingData, testData, trainingLabels, testLabels] = makeSVMmatrices(parName, segPoints, trainSet, testSet, domSegMid_sim, domSegMid_riv, mixSegProp, teStartPt, teEndPt, rivLat, numStates, paramsFile);
            
            % save participant testing and training sets
            save(['ML/parTrainSets/' parName '_' trainSet '_' num2str(numStates) 'states'], 'trainingData', 'trainingLabels');
            if strcmp(testSet, 'identical')
                save(['ML/parTestSets/' parName '_' trainSet '_' num2str(numStates) 'states'], 'testData', 'testLabels');
            else
                save(['ML/parTestSets/' parName '_' testSet '_' num2str(numStates) 'states'], 'testData', 'testLabels');
            end
        end
        
        % Ensure that each class has the same number of training examples
        if withMixed
            if numStates == 2
                % how many of each class are there for this participant?
                mixedRows = find(strcmp(trainingLabels, 'mi'));
                numMixed = length(mixedRows);
                domRows = find(strcmp(trainingLabels, 'dom'));
                numDom = length(domRows);
                
                % randomly choose training examples to keep
                if numDom > numMixed
                    keepInds = randi(numDom, 1, numMixed);
                    keepDomRows = domRows(keepInds);
                    keepRows = [mixedRows; keepDomRows];
                elseif numMixed > numDom
                    keepInds = randi(numMixed, 1, numDom);
                    keepDomRows = mixedRows(keepInds);
                    keepRows = [domRows; keepDomRows];
                else
                    keepRows = [mixedRows; domRows];
                end
                
            elseif numStates == 3
                % how many of each class are there for this participant?
                loRows = find(strcmp(trainingLabels, 'lo'));
                numLo = length(loRows);
                hiRows = find(strcmp(trainingLabels, 'hi'));
                numHi = length(hiRows);
                miRows = find(strcmp(trainingLabels, 'mi'));
                numMi = length(miRows);
                
                % randomly choose training examples to keep
                minNum = min([numLo numHi numMi]);
                if numLo > minNum
                    keepInds = randi(numLo, 1, minNum);
                    keepLoRows = loRows(keepInds);
                else
                    keepLoRows = loRows;
                end
                if numHi > minNum
                    keepInds = randi(numHi, 1, minNum);
                    keepHiRows = hiRows(keepInds);
                else
                    keepHiRows = hiRows;
                end
                if numMi > minNum
                     keepInds = randi(numMi, 1, minNum);
                    keepMiRows = miRows(keepInds);
                else
                    keepMiRows = miRows;
                end
                
                keepRows = [keepLoRows; keepHiRows; keepMiRows];
            end
            
            keepRows = sort(keepRows);
            trainingLabels = trainingLabels(keepRows,:);
            trainingData = trainingData(keepRows,:);
        end
        
        % Collect training and testing data for all participants
        groupTrainingData = [groupTrainingData; trainingData];
        groupTrainingLabels = [groupTrainingLabels; trainingLabels];
        groupTestData = [groupTestData; testData];
        groupTestLabels = [groupTestLabels; testLabels];
        trainParVect = [trainParVect; repmat(parNum, length(trainingLabels), 1)];
        testParVect = [testParVect; repmat(parNum, length(testLabels), 1)];
    end
end

%% Standardize training and testing datasets
groupTrainingData = zscore(groupTrainingData, 0, 1);
troupTestData = zscore(groupTestData, 0, 1);

%% Iterate through potential models

for iModel = 1:length(models)
    model = models{iModel};
    
    if strcmp(model, 'linear')
        modSigmas = NaN;
    else
        modSigmas = sigmas;
    end
    
    %% Iterate through potential kernel scales
    
    tuningMat = nan(length(modSigmas),length(Cs));
    resubMat = nan(length(modSigmas),length(Cs));
    confMat = cell(length(modSigmas),length(Cs));
    
    for iKs = 1:length(modSigmas)
        ksVal = modSigmas(iKs);
        
        %% Iterate through potential box constraints
        
        for iC = 1:length(Cs);
            cVal = Cs(iC);
            
            %% Space to store results
            
            accuracies = [];
            resubAccs = [];
            postProbs = [];
            
            if withMixed
                confusionMats = {};
            end
            %mixedAccuracies = [];
            
            %% Classify dominant percepts for each participant by training on the remaining participants
            
            Msglength = 0;
            
            for iGroup = 1:length(groupParNums)
                group = groupCodes{iGroup};
                for iTestPar = 1:length(groupParNums{iGroup})
                    
                    fprintf(repmat('\b',1, Msglength));
                    Msg = sprintf('Running participant %d of %d\n', (iGroup - 1)*length(groupParNums{1}) + iTestPar, length(groupParNums{1}) + length(groupParNums{2}));
                    fprintf(Msg);
                    Msglength = numel(Msg);
                    
                    % participant-specific set-up
                    testParNum = groupParNums{iGroup}(iTestPar);
                    numFormat = '%02d';
                    testPar = [group num2str(testParNum, numFormat)];
                    
                    % Get EEG file info to determine whether participant exists
                    EEGfiles = dir([eegDir testPar '*']);
                    % ...if the participant exists...
                    if isempty(EEGfiles)
                        continue;
                    end
                    
                    
                    %% Train the classifier
                    
                    trData = groupTrainingData((trainParVect ~= testParNum),:);
                    trLabels = groupTrainingLabels((trainParVect ~= testParNum),:);
                    
                    if withMixed
                        if numStates == 3
                            %t = templateSVM('Standardize', false, 'BoxConstraint', cVal);
                            %SVMModel = fitcecoc(trData, trLabels, 'Learners', t, 'ClassNames', {'lo', 'mi', 'hi'});
                            SVMModel = TreeBagger(100, trData, trLabels);
                            
                            % 10-fold cross-validation
                            %SVMModel = crossval(Mdl);
                        elseif numStates == 2
                            if strcmp(model, 'linear')
                                SVMModel = fitcsvm(trData, trLabels, 'KernelFunction', 'linear', 'BoxConstraint', cVal, 'Standardize', false, 'ClassNames', {'dom','mi'});
                            elseif strcmp(model, 'rbf')
                                SVMModel = fitcsvm(trData, trLabels, 'KernelFunction', 'rbf', 'BoxConstraint', cVal, 'KernelScale', ksVal, 'Standardize', false, 'ClassNames', {'dom','mi'});
                                %SVMModel = fitcsvm(trData, trLabels, 'KernelFunction', 'rbf', 'BoxConstraint', 1, 'Standardize', false, 'ClassNames', {'dom','mi'});
                            end
                        end
                    else
                        SVMModel = fitcsvm(trData, trLabels, 'KernelFunction', 'linear', 'Standardize', true, 'ClassNames', {'lo','hi'});
                    end
                    
                    %% Classify test data
                    
                    teData = groupTestData((trainParVect == testParNum),:);
                    teLabels = groupTestLabels((trainParVect == testParNum),:);
                    [classedLabels, scores] = predict(SVMModel, teData);
                    [resubLabels, ~] = resubPredict(SVMModel);
                    
                    %% Check classifications, store accuracies
                    
                    correct = nan(length(teLabels),1);
                    resubCorrect = nan(length(teLabels),1);
                    %mixedCorrect = [];
                    for i = 1:length(teLabels)
                        if withMixed
                            if numStates == 3
                                if strcmp(teLabels(i), 'hi') || strcmp(teLabels(i), 'lo') || strcmp(teLabels(i), 'mi')
                                    correct(i) = strcmp(teLabels(i), classedLabels(i));
                                    resubCorrect(i) = strcmp(trLabels(i), resubLabels(i));
                                end
                                
                                % confusion matrix
                                C = confusionmat(teLabels, classedLabels, 'order', {'lo','mi','hi'});
                                
                                %             % mixed state accuracy
                                %             if strcmp(teLabels(i), 'mi')
                                %                 mixedCorrect = [mixedCorrect correct(i)];
                                %             end
                            elseif numStates == 2
                                if strcmp(teLabels(i), 'dom') || strcmp(teLabels(i), 'mi')
                                    correct(i) = strcmp(teLabels(i), classedLabels(i));
                                    resubCorrect(i) = strcmp(trLabels(i), resubLabels(i));
                                end
                                
                                % confusion matrix
                                C = confusionmat(teLabels, classedLabels, 'order', {'dom','mi'});
                            end
                        else
                            % ignore mixed states
                            if strcmp(teLabels(i), 'hi') || strcmp(teLabels(i), 'lo')
                                correct(i) = strcmp(teLabels(i), classedLabels(i));
                            end
                        end
                    end
                    
                    scores = max(scores,[], 2);
                    
                    accuracies = [accuracies; nanmean(correct)];
                    resubAccs = [resubAccs; nanmean(resubCorrect)];
                     postProbs = [postProbs; nanmean(scores)];
                    
                    if withMixed
                        confusionMats = [confusionMats C];
                    else
                        confusionMats = [];
                    end
                    %mixedAccuracies = [mixedAccuracies; nanmean(mixedCorrect)];
                end
            end
            tuningMat(iKs, iC) = nanmean(accuracies);
            resubMat(iKs, iC) = nanmean(resubAccs);
            confMat{iKs, iC} = confusionMats;
        end
    end
    
    %% Plotting
    
    if strcmp(model, 'linear')
        figure
        subplot(1,2,1)
        plot(Cs, tuningMat);
        set(gca,'XScale','log');
        xlabel('Box Constraint')
        title('Cross-Validation Accuracies')
        
        subplot(1,2,2)
        plot(Cs, resubMat);
        set(gca, 'XScale', 'log');
        xlabel('BoxConstraint')
        title('Resubstitution Accuracies')
        
    elseif strcmp(model, 'rbf')
        figure
        subplot(1,2,1)
        surf(sigmas, Cs, zeros(size(tuningMat)), tuningMat);
        view(0,90);
        set(gca,'XScale','log');
        set(gca, 'YScale', 'log');
        xlabel('Box Constraint');
        ylabel('Sigma');
        colorbar
        title('Cross-Validation Accuracies')

        subplot(1,2,2)
        surf(sigmas, Cs, zeros(size(tuningMat)), resubMat);
        view(0,90);
        set(gca,'XScale','log');
        set(gca, 'YScale', 'log');
        xlabel('Box Constraint');
        ylabel('Sigma');
        colorbar
        title('Resubstitution Accuracies')
    end
end
%end

%toc
%end