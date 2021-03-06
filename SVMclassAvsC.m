% Uses SVM to classify participants' diagnosis
% uses Matlab's SVM toolbox

%loadPaths % add dependency folders to path

clearvars
tic
maxResamples = 1; %Ideally ~10,000 but try with 10 or 100
%If you run this 10000, you'll get 70% accuracy in diagnosing
%tic

% load parameters
%run(paramsFile);
%analysisParams
%groupAnalysisParams

%% Parameters
% Location of Raw Data
if exist('../../runScripts/rawData', 'dir')
    eegDir = '../../runScripts/rawData/EEG/'; 
    keyPressDir = '../../runScripts/rawData/Behavior/'; 
else
    error('please provide directory where raw data is stored');
end

groupParNums = {[105:138], [1:28]}; %{[105:138], [1:28]}; %{[106 112 131 132], []}; %{[105:124 126:138], [1:29]};%{[105:138], [01:13 15:20]}; %, 01:20};
groupCodes = {'stratus','cumulus'}; %{'stratus','cumulus'};

% exclude participants
load('pre-processing/stratExcluded.mat')
load('pre-processing/cumExcluded.mat')
groupParNums{1} = setdiff(groupParNums{1},stratExcluded);
groupParNums{2} = setdiff(groupParNums{2},cumExcluded);

runTypes = {'dartRival','dartSim'}; %'dartSim', 'marzRival'
runIndices = [1:3]; % generally 1:3

stimFreqs = [5.666 8.5]; % hz
harFreqs = [2*stimFreqs(1) 2*stimFreqs(2)]; % harmonics
imFreqs = [2*stimFreqs(2) - stimFreqs(1) stimFreqs(1) + stimFreqs(2)]; % intermodulation frequencies

%% 
predictors = [];
parVect = [];
diagnoses = {};

%% Load predictors for all participants

for iGroup = 1:length(groupParNums)
    group = groupCodes{iGroup};
    for parNum = groupParNums{iGroup}
     
                
        %% Participant-specific set-up
        
        numFormat = '%02d';
        parName = [group num2str(parNum, numFormat)];
        
        % Get EEG file info
        EEGfiles = dir([eegDir parName '*']);
        % ...if the participant exists...
        if isempty(EEGfiles)
            continue;
        end
        
        parPredictors = [];
        
        
        %% FFTs of RLS Data
        %load(['OscFFTs/' parName '.mat']);
        %predictors = [parPredictors meanfft(freq <= 0.3)];
        
        %% Get OscFreqs
        load(['indices/rateFreqs/' parName '.mat']); 
        parPredictors = [parPredictors oscFreq];
        
        %% Relative FFT amplitudes
        
        electrodes = [29];
        FOIs = [stimFreqs harFreqs imFreqs]; %[stimFreqs harFreqs imFreqs];
        %FOIs = [stimFreqs]; %[stimFreqs harFreqs imFreqs];
        for elec = electrodes
            for iRunType = 1 % dartRival only
                runType = runTypes{iRunType};
                parFFTs = [];
                for runIndex = runIndices
                    runName = [runType num2str(runIndex)];
                    if exist(['FFTs/' parName '_' runName '_' num2str(elec) '.mat'], 'file') % won't exist if run was missing, or all trials in run had low SNR
                        load(['FFTs/' parName '_' runName '_' num2str(elec) '.mat'])
                        if size(elecFFT, 2) == 1 % only 1 good trial in this run
                            elecFFT = elecFFT';
                            parFFTs = [parFFTs; (elecFFT)];
                        else
                        parFFTs = [parFFTs; (elecFFT)]; %aggregates across trials
                        end
                    end
                end
                
            end
            
            % average FFT across trials
            avgFFT = nanmean(parFFTs, 1); 

            % Get frequencies of interest
            FOImat = repmat(FOIs, length(freqAxis), 1);
            freqAxMat = repmat(freqAxis', 1, length(FOIs));
            [~, FOIinds] = min(abs(freqAxMat - FOImat));
                 
            %Make predictor
            %FOIfeat = avgFFT(FOIinds)
                        
            %FOIfeat = avgFFT(FOIinds) - nanmean(avgFFT(FOIinds));
            FOIfeat = zscore(avgFFT(FOIinds)); %warning: you can't zscore only two numbers (imFreqs)
            %FOIfeat = nanmean(avgFFT(FOIinds(3:4))) / nanmean(avgFFT(FOIinds(1:2)));
            parPredictors = [parPredictors FOIfeat];   
            
            clear FOIfeat elecFFT
        end
        
        predictors = [predictors; parPredictors];
        
        % record participants and diagnosis
        parVect = [parVect parNum];
        if strcmp(group, 'stratus')
            diagnoses = [diagnoses; 'C'];
        elseif strcmp(group, 'cumulus')
            diagnoses = [diagnoses; 'A'];
        else
            error('Make sure groups include only stratus and cumulus');
        end
        
        clear parPredictors
    end
    
    
end

%% Standardize training and testing datasets

% z-score each column (feature)
%predictors(2:end,:) = zscore(predictors(2:end,:), 0, 1);

%% Train the classifier

accuracies = [];
postProbs = [];
Msglength = 0;

firstDiagnoses = diagnoses;
firstPredictors = predictors;
firstParVect = parVect;

for iResample = 1:maxResamples
    
    groupNums = find((strcmp(firstDiagnoses,'C')));
    nonGroupNums = find((strcmp(firstDiagnoses,'A')));
    if maxResamples > 1
        %Randomly delete one control to even group #
        out = randi(length(groupNums));
        randIn = setdiff(groupNums,out);
        randIn = [randIn; nonGroupNums];

        diagnoses = firstDiagnoses(randIn,:); 
        predictors = firstPredictors(randIn,:);
        parVect = firstParVect(:,randIn);
    end
    
    keepDiagnoses = diagnoses;
    keepPredictors = predictors;
    keepParVect = parVect;
    
    
    %Msg = sprintf('Running iteration %d of %d\n', iResample, maxResamples);
    %fprintf(Msg);
    
    [sum(strcmp(diagnoses,'C')) sum(strcmp(diagnoses,'A'))];

    keepclassedLabel = [];
    for iTestPar = 1:(length(keepDiagnoses)-1)
        
        diagnoses = keepDiagnoses;
        predictors = keepPredictors;
        parVect = keepParVect;

        fprintf(repmat('\b',1, Msglength));
        %Msg = sprintf('Running participant %d of %d\n', (iGroup - 1)*length(groupParNums{1}) + iTestPar, length(groupParNums{1}) + length(groupParNums{2}));
        %fprintf(Msg);
        %Msglength = numel(Msg);

        % participant-specific set-up
        testParNum = parVect(iTestPar);
        if strcmp(diagnoses(iTestPar),'C')
            iGroup = 1;
        else strcmp(diagnoses(iTestPar),'A')
            iGroup = 2;
        end

        group = groupCodes{iGroup};

        numFormat = '%02d';
        testPar = [group num2str(testParNum, numFormat)];

        % Get EEG file info to determine whether participant exists
        EEGfiles = dir([eegDir testPar '*']);
        % ...if the participant exists...
        if isempty(EEGfiles)
            continue;
        end

        % CER 2017 - is it OK that this happens?
        if isempty(intersect(parVect,testParNum))
            hhhh
            continue
        end


        %% Train the classifier
        trData = predictors((parVect ~= testParNum),:);
        trLabels = diagnoses((parVect ~= testParNum),:);
        diagnoses = diagnoses((parVect ~= testParNum),:); %% Cer 2017 to resample               
        parVect = parVect(parVect ~= testParNum);
        
        %[sum(strcmp(trLabels,'C')) sum(strcmp(trLabels,'A'))]

        if maxResamples > 1
            %% Cer 2017 randomly removes one subject to even groups  
            if iGroup == 2
                groupNums = find((strcmp(diagnoses,'C')));
                nonGroupNums = find((strcmp(diagnoses,'A')));
                randOut = randi(length(groupNums));
                randIn = setdiff(groupNums,randOut);
                randIn = [randIn; nonGroupNums];
            elseif iGroup == 1
                groupNums = find((strcmp(diagnoses,'A')));
                nonGroupNums = find((strcmp(diagnoses,'C')));
                randOut = length(nonGroupNums) + randi(length(groupNums));
                randIn = setdiff(groupNums,randOut);
                randIn = [nonGroupNums; randIn];
            end

            trData = trData(randIn,:);
            trLabels = trLabels(randIn,:);
            diagnoses = diagnoses(randIn,:);
            predictors = predictors(randIn,:);
            parVect = parVect(:,randIn);

            [sum(strcmp(trLabels,'C')) sum(strcmp(trLabels,'A'))];

            if sum(strcmp(trLabels,'C')) ~= sum(strcmp(trLabels,'A'))
                hhh %Groups need to be matched for size!
            end
            
        
            %Resample data (bootstrap)
            nGroup = sum(strcmp(trLabels,'C'));
            cPerm = randi(nGroup,nGroup,1);
            aPerm = randi(nGroup,nGroup,1)+nGroup;
            permI = [cPerm; aPerm];

            diagnoses = diagnoses(permI,:); 
            predictors = predictors(permI,:);
            parVect = parVect(:,permI);
            trData = trData(permI,:);
            trLabels = trLabels(permI,:);
            
            
        end
        
        [sum(strcmp(trLabels,'C')) sum(strcmp(trLabels,'A'))];
        SVMModel = fitcsvm(trData, trLabels,  ... 
            'Standardize', true, 'ClassNames', {'C','A'}, ...
            'KernelFunction','linear');
        ScoreSVMModel = fitSVMPosterior(SVMModel);
        
%         sv = SVMModel.SupportVectors;
%         figure
%         gscatter(trData(:,1),trData(:,2),trLabels)
%         hold on
%         plot(sv(:,1),sv(:,2),'ko','MarkerSize',10)
%         legend('versicolor','virginica','Support Vector')
%         hold off
        
        

        %% Classify test data
        teData = keepPredictors((keepParVect == testParNum),:);
        teLabel = keepDiagnoses((keepParVect == testParNum),:);
        [classedLabel, score] = predict(ScoreSVMModel, teData);
        resubLabels = resubPredict(SVMModel);

        %% Check classifications, store accuracies and posterior probabilities
        correct = strcmp(teLabel, classedLabel);        
        score = max(score,[], 2);

        accuracies = [accuracies; correct];
        postProbs = [postProbs; score];
        keepclassedLabel = [keepclassedLabel classedLabel];

        clear resubLabels teLabel teData trData trLabels correct score
        
    end
    keepclassedLabel;
    
    meanAcc(iResample) = nanmean(accuracies);
    [iResample]
end

%Print to screen
disp(['mean: ' num2str(mean(meanAcc))])
disp(['std: ' num2str(std(meanAcc))])

%% Plot and Save Figure
accMean = nanmean(meanAcc);
accError = std(meanAcc);

setFigProps
figure
hold on
bar(1,[accMean])
Labels = {'Accuracy'};
set(gca, 'XTick', 1, 'XTickLabel', Labels);
errorbar(1,[accMean  ],[accError ],'.')
axis([0 2 0.25 .9])
hline = refline(0,0.5);
hline.Color = 'r';
title(['Classify Diagnosis, Combined Groups  .'])
ylabel('Percent Correct %')

%saveas(gcf,'plots/barGraphs/Classify_Diagnosis','jpg')


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

toc