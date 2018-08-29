% Train an SVM to classify dominant state of a left out participant using
% button presses of all other participants 
% Dependencies: svmClassLeftOutPar

loadPaths % add dependency folders to path

% clearvars

tic

% load parameters (includes SVM features)
groupAnalysisParams
svmParams

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
    
    if acrossPars
        [accuracies, postProbs, confusionMats, mixedAccs, domAccs] = svmClassLeftOutPar(trainSet, testSet, rivLat, segPoints, domSegMid_sim, domSegMid_riv, mixSegProp, teStartPt, teEndPt, withMixed, numStates);
    elseif withMixed
        [accuracies, postProbs, confusionMats, mixedAccs, domAccs] = svmClassify_withMixed(trainSet, testSet, rivLat, segPoints, domSegMid_sim, domSegMid_riv, mixSegProp, teStartPt, teEndPt, numStates);
    else
        [accuracies, postProbs] = svmClassify(trainSet, testSet, rivLat, segPoints, domSegMid_sim, domSegMid_riv, mixSegProp, teStartPt, teEndPt, numStates);
        confusionMats = NaN;
        mixedAccs = NaN;
        domAccs = NaN;
    end
    
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
        if exist('mixedAccs', 'var')
            disp(['MIXED mean: ' num2str(nanmean(mixedAccs))])
        end
        if exist('domAccs', 'var')
            disp(['DOM mean: ' num2str(nanmean(domAccs))])
        end
        disp(['']);
        
        figure
        histogram(accuracies, 15);
        xlabel('accuracy')
        if acrossPars
            title('Accuracies: classify a left out participant       .')
        else
            title('Accuracies: classify within participant      .')
        end
    end
end

% save accuracies & posterior probabilities timecourse for each participant
if acrossPars
    saveLabel = 'acrossPars';
else
    saveLabel = 'withinPars';
end
if ~isempty(strfind(testSet, 'Transitions'))
    save([saveDir 'run' runID '_' group '_train-' trainSet '_test-' testSet '_' num2str(length(accuracies)) 'pars_' classes '_' saveLabel], 'accuracyMat', 'postProbMat', 'domAccMat', 'mixAccMat', 'teTime');
end
fprintf('Saving accuracies & probabilities timecourse for each participant\n')

toc
