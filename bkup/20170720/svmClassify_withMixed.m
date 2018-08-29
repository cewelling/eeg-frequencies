function [accuracies, postProbs, confusionMats, mixedAccs, domAccs] = svmClassify_withMixed(trainSet, testSet, rivLat, segPoints, domSegMid_sim, domSegMid_riv, mixSegProp, teStartPt, teEndPt, numStates, paramsFile)

% Uses classifier to...
% uses Matlab's SVM toolbox

%loadPaths % add dependency folders to path

%clearvars

%tic

% load parameters
run(paramsFile);

accuracies = [];
postProbs = [];
confusionMats = {};

mixedAccs = [];
domAccs = [];

%% Iterate through participants

for iPar = 1:length(parNums)
    
    %% Participant-specific set-up
    
    numFormat = '%02d';
    parName = [group num2str(parNums(iPar), numFormat)];
    
    % Get EEG file info
    EEGfiles = dir([eegDir parName '*']);
    % ...if the participant exists...
    if isempty(EEGfiles)
        continue;
    end
    
    %% Create matrices of training and test data
    
    [trainingData, testData, trainingLabels, testLabels, labelVect] = makeSVMmatrices(parName, segPoints, trainSet, testSet, domSegMid_sim, domSegMid_riv, mixSegProp, teStartPt, teEndPt, rivLat, numStates, paramsFile);
    
    % find NaN rows 
    trainingData_all = trainingData;
    trainingData = trainingData_all(~isnan(nanmean(trainingData_all, 2)),:);
    trainingLabels = trainingLabels(~isnan(nanmean(trainingData_all, 2)),:);
    testData_all = testData;
    testData = testData_all(~isnan(nanmean(testData_all, 2)),:);
    testLabels = testLabels(~isnan(nanmean(testData_all,2)),:);
    
    m = 1; % max ratio of examples in one class to another
    
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
    
    % SNR not high enough; skip this participant
    if isempty(trainingData) || isempty(testData)
        continue;
    end
    
    % Not enough of one class
    if size(trainingData, 1) < 30
        continue;
    end
    
    %% Train the classifier
    
    t = templateSVM('Standardize', 1, 'BoxConstraint', 0.01);
    
    if numStates == 3
        SVMModel = fitcecoc(trainingData, trainingLabels, 'Learners', t, 'fitPosterior', 1, 'ClassNames', {'lo', 'mi', 'hi'});
    elseif numStates == 2
        SVMModel = fitcsvm(trainingData, trainingLabels, 'KernelFunction', 'linear', 'Standardize', true, 'ClassNames', {'dom','mi'});
    end
    
    %SVMModel = fitcsvm(trainingData, trainingLabels, 'KernelFunction', 'linear', 'Standardize', true, 'ClassNames', {'lo','hi'});
    %SVMModel = fitcsvm(trainingData, trainingLabels, 'KernelFunction', 'rbf', 'BoxConstraint', 1, 'Standardize', true, 'ClassNames', {'lo','hi'});
    %SVMModel = fitcsvm(trainingData, trainingLabels);
    %%
    %         %% Visualize the training dataset (if we only have 2 dimensions)
    %
    %             sv = SVMModel.SupportVectors;
    %             figure
    %             gscatter(trainingData(:,1),trainingData(:,2),trainingLabels,[],[],[],'on')
    %             hold on
    %             %plot(sv(:,1),sv(:,2),'ko','MarkerSize',10) % only works if we don't standardize before training
    %             %legend('low dominant','high dominant') %,'support vector')
    %             title(['Training on ' trainSet ': ' parName '    .']);
    %             hold off
    %
    %         %% Visualize the testing dataset (if we only have 2 dimensions)
    %
    %         sv = SVMModel.SupportVectors;
    %         figure
    %         hiOrLo = find(strcmp(testLabels, 'hi') | strcmp(testLabels, 'lo'));
    %         gscatter(testData(hiOrLo,1),testData(hiOrLo,2),testLabels(hiOrLo),[],[],[],'on')
    %         hold on
    %         %plot(sv(:,1),sv(:,2),'ko','MarkerSize',10) % only works if we don't standardize before training
    %         %legend('low dominant','high dominant') %,'support vector')
    %         title(['Testing on ' testSet ': ' parName '    .']);
    %         hold off
    
    %% Classify test data
    
    [classedLabels, scores] = predict(SVMModel, testData);
    
    %% Check classifications, store accuracies and posterior probabilities
    
    correct = nan(length(testLabels),1);
    mixedCorrect = [];
    domCorrect = [];
    
    for i = 1:length(testLabels)
        
        %             if ~strcmp(testLabels(i), 'mi') && ~isempty(testLabels{i}) % ignore mixed states for the moment, bins with multiple states
        %                 correct(i) = strcmp(testLabels(i), classedLabels(i));
        %             end
        
        if numStates == 3
            if strcmp(testLabels(i), 'hi') || strcmp(testLabels(i), 'lo') || strcmp(testLabels(i), 'mi')
                correct(i) = strcmp(testLabels(i), classedLabels(i));
            end
            
            % separate accuracies for each state
            if strcmp(testLabels(i), 'mi')
                mixedCorrect = [mixedCorrect correct(i)];
            elseif strcmp(testLabels(i), 'hi') || strcmp(testLabels(i), 'lo')
                domCorrect = [domCorrect correct(i)];
            end
            
        elseif numStates == 2
            if strcmp(testLabels(i), 'dom') || strcmp(testLabels(i), 'mi')
                correct(i) = strcmp(testLabels(i), classedLabels(i));
            end
        end
    end
    
    scores = max(scores,[], 2);
    
    % confusion matrix
    if numStates == 3
        C = confusionmat(testLabels, classedLabels, 'order', {'lo','mi','hi'});
    elseif numStates == 2
        C = confusionmat(testLabels, classedLabels, 'order', {'dom','mi'});
    end
    
    accuracies = [accuracies; nanmean(correct)];
    postProbs = [postProbs; nanmean(scores)];
    confusionMats = [confusionMats; C];
    
    mixedAccs = [mixedAccs; nanmean(mixedCorrect)];
    domAccs = [domAccs; nanmean(domCorrect)];
    
end
end

%     % save accuracies
%     save(['ML/segAccs/' group '_train-' trainSet '_test-' testSet '_' num2str(size(testData, 2)) 'features'], 'accuracies')
%
%     disp(['group: ' group])
%     %disp(['time after epoch start: ' num2str(segMiddle) ' s'])
%     disp(['mean: ' num2str(nanmean(accuracies))])
%     disp(['median: ' num2str(nanmedian(accuracies))])
%     disp(['min: ' num2str(nanmin(accuracies))])
%     disp(['max: ' num2str(nanmax(accuracies))])
%     disp(['']);
%
%     figure
%     histogram(accuracies, 15);
%     xlabel('accuracy')

%toc
%end