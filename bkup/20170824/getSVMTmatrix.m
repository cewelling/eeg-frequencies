function [ lowMat, highMat, labelMat ] = getTmatrix( parName, date, runType, channel, freqType, rivLat)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% load parameters
groupAnalysisParams;
svmParams

zScore = 0; % z-score transitions?

lowMat = [];
highMat = [];
labelMat = [];

rivLat_str = num2str(round(rivLat*1000)); % for use in filenames

for iRunType = runType % 1 for dartRival, 2 for dartSim
    cRunType = runTypes{iRunType};
    for runIndex = runIndices
        runName = [cRunType num2str(runIndex)];
        
        % Get the appropriate EEG file for this run
        EEGfile = [eegDir parName '_' runName '_' date '.bdf'];
        
        % skip over participants / runs that don't exist
        if ~exist(EEGfile, 'file')
            continue;
        end
        
        % load participant's transitions
        if strcmp(cRunType, 'dartSim') && useSimSchedule
            if ~exist(['transitions/transitionsByRun_simSched/' parName '_' runName '_' normType '_' freqType{1} '_' num2str(channel) '.mat'], 'file')
                makeParPlot(parName, date);
            end
            load(['transitions/transitionsByRun_simSched/' parName '_' runName '_' normType '_' freqType{1} '_' num2str(channel) '.mat'])
        else
            if ~exist(['transitions/transitionsByRun/' parName '_' runName '_' normType '_' freqType{1} '_' num2str(channel) '.mat'], 'file')
                makeParPlot(parName, date);
            end
            load(['transitions/transitionsByRun/' parName '_' runName '_' normType '_' freqType{1} '_' num2str(channel) '.mat'])
        end
        
        % matrices of low to high transitions
        l2hLowMat_run = run_tLists{1}.f1;
        l2hHighMat_run = run_tLists{1}.f2;
        % label matrices
        l2hLabelMat_run = getSegPresses(parName, runName, run_tLists{1}.tTrials, run_tLists{1}.tPoints, rivLat);
        
        % matrices of high to low transitions
        h2lLowMat_run = run_tLists{2}.f1;
        h2lHighMat_run = run_tLists{2}.f2;
        % label matrices
        h2lLabelMat_run = getSegPresses(parName, runName, run_tLists{2}.tTrials, run_tLists{2}.tPoints, rivLat);
        
        % Combine matrices across runs and transition types
        lowMat = [lowMat; l2hLowMat_run; h2lLowMat_run];
        highMat = [highMat; l2hHighMat_run; h2lHighMat_run];
        labelMat = [labelMat; l2hLabelMat_run; h2lLabelMat_run];
    end
end

% Z-score transitions
if zScore == 1
    lowMat = zscore(lowMat, 0, 2);
    highMat = zscore(highMat, 0, 2);
end

% save matrices of transitions and corresponding labels
if runType == 1 % rivalry transition
    save(['ML/tMatrices/' parName '_' freqType{1} '_' num2str(channel) '_' runTypes{runType} '_' rivLat_str 'msLat'], 'lowMat', 'highMat', 'labelMat');
else %sim transition
    save(['ML/tMatrices/' parName '_' freqType{1} '_' num2str(channel) '_' runTypes{runType}], 'lowMat', 'highMat', 'labelMat');
end
fprintf(['saving transition matrices and corresponding labels for ' parName '\n']);

end

