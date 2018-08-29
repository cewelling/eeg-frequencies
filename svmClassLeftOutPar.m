function [accuracies, postProbs, confusionMats, mixedAccs, domAccs] = svmClassLeftOutPar(trainSet, testSet, rivLat, segPoints, domSegMid_sim, domSegMid_riv, mixSegProp, teStartPt, teEndPt, withMixed, numStates)

% Uses classifier to...
% uses Matlab's SVM toolbox

%loadPaths % add dependency folders to path

%clearvars

%tic

% load parameters
groupAnalysisParams
svmParams

accuracies = [];
postProbs = [];

if withMixed
    confusionMats = {};
end

mixedAccs = [];
domAccs = [];

%% Iterate through participants to create training and testing sets

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
        elseif strfind(testSet, 'Transitions')
            testExist = exist(['ML/parTestSets/' parName '_' testSet '_' num2str(numStates) 'states_' num2str(teStartPt) 'to' num2str(teEndPt) '.mat'], 'file');
            if testExist
                load(['ML/parTestSets/' parName '_' testSet '_' num2str(numStates) 'states_' num2str(teStartPt) 'to' num2str(teEndPt) '.mat']);
            end
        end
        
        if ~(trainExist && testExist)
            % Get training and testing data for this participant
            [trainingData, testData, trainingLabels, testLabels] = makeSVMmatrices(parName, segPoints, trainSet, testSet, domSegMid_sim, domSegMid_riv, mixSegProp, teStartPt, teEndPt, rivLat, numStates);
            
            % save participant testing and training sets
            save(['ML/parTrainSets/' parName '_' trainSet '_' num2str(numStates) 'states'], 'trainingData', 'trainingLabels');
            if strcmp(testSet, 'identical')
                save(['ML/parTestSets/' parName '_' trainSet '_' num2str(numStates) 'states'], 'testData', 'testLabels');
            elseif strfind(testSet, 'Transitions')
                save(['ML/parTestSets/' parName '_' testSet '_' num2str(numStates) 'states_' num2str(teStartPt) 'to' num2str(teEndPt)], 'testData', 'testLabels');
            end
            fprintf('Saving participant testing and training sets\n');
        end
        
        % find NaN rows (transitions that were cut off)
        trainingData_all = trainingData;
        trainingData = trainingData_all(~isnan(nanmean(trainingData_all, 2)),:);
        trainingLabels = trainingLabels(~isnan(nanmean(trainingData_all, 2)),:);
        testData_all = testData;
        testData = testData_all(~isnan(nanmean(testData_all, 2)),:);
        testLabels = testLabels(~isnan(nanmean(testData_all,2)),:);
        
        m = 1; % max ratio of examples in one class to another;
        % so if this is 1, the two classes must have an equal # of training examples
        
        % Ensure that each class has the same number of training examples
        if withMixed
            if numStates == 2
                % how many of each class are there for this participant?
                mixedRows = find(strcmp(trainingLabels, 'mi'));
                numMixed = length(mixedRows);
                domRows = find(strcmp(trainingLabels, 'dom'));
                numDom = length(domRows);
                
                % randomly choose training examples to keep
                if numDom > numMixed*m
                    keepInds = randi(numDom, 1, numMixed*m);
                    keepDomRows = domRows(keepInds);
                    keepRows = [mixedRows; keepDomRows];
                elseif numMixed > numDom*m
                    keepInds = randi(numMixed, 1, numDom*m);
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
                minNum = min([numLo numHi numMi])*m;
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

% % find NaN rows (transitions that were cut off)
% groupTrData_all = groupTrainingData;
% groupTrainingData = groupTrData_all(~isnan(nanmean(groupTrData_all, 2)),:);
% groupTrainingLabels = groupTrainingLabels(~isnan(nanmean(groupTrData_all, 2)),:);
% trainParVect = trainParVect(~isnan(nanmean(groupTrData_all, 2)),:);
% groupTeData_all = groupTestData;
% groupTestData = groupTeData_all(~isnan(nanmean(groupTeData_all, 2)),:);
% groupTestLabels = groupTestLabels(~isnan(nanmean(groupTeData_all,2)),:);
% testParVect = testParVect(~isnan(nanmean(groupTeData_all, 2)),:);

% z-score each column (feature)
groupTrainingData = zscore(groupTrainingData, 0, 1);
groupTestData = zscore(groupTestData, 0, 1);

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
                t = templateSVM('Standardize', 1, 'BoxConstraint', boxConstraint);
                SVMModel = fitcecoc(trData, trLabels, 'Learners', t, 'ClassNames', {'lo', 'mi', 'hi'});
                %SVMModel = TreeBagger(100, trData, trLabels);
                
                if plot2features
                    % Plot 2 features
                    figure
                    hiOrLo = find(strcmp(trLabels, 'hi') | strcmp(trLabels, 'lo') | strcmp(trLabels, 'mi'));
                    gscatter(trData(hiOrLo,75),trData(hiOrLo,76),trLabels(hiOrLo),[],[],[],'on')
                    hold on
                    %plot(sv(:,1),sv(:,2),'ko','MarkerSize',10) % only works if we don't standardize before training
                    %legend('low dominant','high dominant') %,'support vector')
                    title(['Training on ' trainSet ': ' testPar '    .']);
                    hold off
                end
                
                % 10-fold cross-validation
                %SVMModel = crossval(Mdl);
            elseif numStates == 2
                SVMModel = fitcsvm(trData, trLabels, 'KernelFunction', 'linear', 'BoxConstraint', boxConstraint, 'Standardize', false, 'ClassNames', {'dom','mi'});
                %SVMModel = fitcsvm(trData, trLabels, 'KernelFunction', 'rbf', 'BoxConstraint', 0.01, 'KernelScale', 1, 'Standardize', false, 'ClassNames', {'dom','mi'});
            end
        else
            SVMModel = fitcsvm(trData, trLabels, 'KernelFunction', 'linear', 'BoxConstraint', boxConstraint, 'Standardize', false, 'ClassNames', {'lo','hi'});
        end
        
        %% Classify test data
        
        teData = groupTestData((testParVect == testParNum),:);
        teLabels = groupTestLabels((testParVect == testParNum),:);
        [classedLabels, scores] = predict(SVMModel, teData);
        
        % Plot first 2 features
        if plot2features
            figure
            hiOrLo = find(strcmp(teLabels, 'hi') | strcmp(teLabels, 'lo') | strcmp(teLabels, 'mi'));
            gscatter(teData(hiOrLo,75),teData(hiOrLo,76),teLabels(hiOrLo),[],[],[],'on')
            hold on
            %plot(sv(:,1),sv(:,2),'ko','MarkerSize',10) % only works if we don't standardize before training
            %legend('low dominant','high dominant') %,'support vector')
            title(['Testing on ' testSet ': ' testPar '    .']);
            hold off
        end
        
        %% Check classifications, store accuracies
        
        correct = nan(length(teLabels),1);
        mixedCorrect = [];
         domCorrect = [];
         for i = 1:length(teLabels)
             if withMixed
                 if numStates == 3
                     if strcmp(teLabels(i), 'hi') || strcmp(teLabels(i), 'lo') || strcmp(teLabels(i), 'mi')
                         correct(i) = strcmp(teLabels(i), classedLabels(i));
                     end
                     
                     % confusion matrix
                     C = confusionmat(teLabels, classedLabels, 'order', {'lo','mi','hi'});
                     
                     % separate accuracies for each state
                     if strcmp(teLabels(i), 'mi')
                         mixedCorrect = [mixedCorrect correct(i)];
                     elseif strcmp(teLabels(i), 'hi') || strcmp(teLabels(i), 'lo')
                         domCorrect = [domCorrect correct(i)];
                     end
                 elseif numStates == 2
                     if strcmp(teLabels(i), 'dom') || strcmp(teLabels(i), 'mi')
                         correct(i) = strcmp(teLabels(i), classedLabels(i));
                     end
                     
                     % confusion matrix
                     C = confusionmat(teLabels, classedLabels, 'order', {'dom','mi'});
                     
                     % separate accuracies for each state
                     if strcmp(teLabels(i), 'mi')
                         mixedCorrect = [mixedCorrect correct(i)];
                     elseif strcmp(teLabels(i), 'dom')
                         domCorrect = [domCorrect correct(i)];
                     end
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
         postProbs = [postProbs; nanmean(scores)];
         
         if withMixed
             confusionMats = [confusionMats C];
         else
             confusionMats = [];
         end
         
         mixedAccs = [mixedAccs; nanmean(mixedCorrect)]; 
         domAccs = [domAccs; nanmean(domCorrect)];
    end
end
end

%toc
%end