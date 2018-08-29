function [accuracies, postProbs] = svmClassify(trainSet, testSet, rivLat, segPoints, domSegMid_sim, domSegMid_riv, mixSegProp, teStartPt, teEndPt, numStates, paramsFile)

% Uses classifier to...
% uses Matlab's SVM toolbox

%loadPaths % add dependency folders to path

%clearvars

%tic

% load parameters
run(paramsFile);

accuracies = [];
postProbs = [];

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
    
    [trainingData, testData, trainingLabels, testLabels] = makeSVMmatrices(parName, segPoints, trainSet, testSet, domSegMid_sim, domSegMid_riv, mixSegProp, teStartPt, teEndPt, rivLat, numStates, paramsFile);
    
    % SNR not high enough; skip this participant
    if isempty(trainingData) || isempty(testData)
        continue;
    end
    
    %% Train the classifier
    
    if numStates == 3
        SVMModel = fitcsvm(trainingData, trainingLabels, 'KernelFunction', 'linear', 'Standardize', true, 'ClassNames', {'lo','hi'});
    elseif numStates == 2
        SVMModel = fitcsvm(trData, trLabels, 'KernelFunction', 'linear', 'Standardize', true, 'ClassNames', {'dom','mi'});
    end
    ScoreSVMModel = fitSVMPosterior(SVMModel);
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
    
    [classedLabels, scores] = predict(ScoreSVMModel, testData);
    
    %% Check classifications, store accuracies and posterior probabilities
    
    correct = nan(length(testLabels),1);
    for i = 1:length(testLabels)
        
        %             if ~strcmp(testLabels(i), 'mi') && ~isempty(testLabels{i}) % ignore mixed states for the moment, bins with multiple states
        %                 correct(i) = strcmp(testLabels(i), classedLabels(i));
        %             end
        
        % ignore mixed states for the moment
        if numStates == 3
            if strcmp(testLabels(i), 'hi') || strcmp(testLabels(i), 'lo')
                correct(i) = strcmp(testLabels(i), classedLabels(i));
            end
        elseif numStates == 2
            if strcmp(testLabels(i), 'dom') || strcmp(testLabels(i), 'mi')
                correct(i) = strcmp(testLabels(i), classedLabels(i));
            end
        end
    end
    
    scores = max(scores,[], 2);
    
    accuracies = [accuracies; nanmean(correct)];
    postProbs = [postProbs; nanmean(scores)];
    
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