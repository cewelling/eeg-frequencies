
clearvars

% load parameters
groupAnalysisParams

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the following matrix for index analysis:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Column 1: Suppression Indices
% Column 2: Slow epoch durations
% Column 3: Fast epoch durations
% Column 4: Mixed epoch durations
% Column 5: Rivalry modulation index (low freq)
% Column 6: Rivalry modulation index (high freq)
% Column 7: Transition rates (transitions per sec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
analMatrix = [];

for iGroup = 1:length(groupParNums)
    for parNum = groupParNums{iGroup}
        numFormat = '%02d';
        parName = [groupCodes{iGroup} num2str(parNum, numFormat)];
        analRow = nan(1, 6);
          
        % Find run indices for runs with desired rivalry stimuli
        if exist([legendDir parName '.txt'],'file')
            Legend = importdata([legendDir parName '.txt']);
        else
            continue;
        end
        runNames = Legend.textdata(:,2);
        searchNames = {};
        for stimNum = stimType
            runName = ['rival' num2str(stimNum)];
            searchNames = [searchNames runName];
        end
        runIndices = [];
        for runIndex = 1:length(runNames)
            for iName = 1:length(searchNames)
                searchName = searchNames{iName};
                if strfind(runNames{runIndex,1}, searchName)
                    runIndices = [runIndices runIndex];
                    break;
                end
            end
        end
        
        %% get suppression indices for rivalry-----------------------------
        
        % get participant files
        parFiles = dir([indicesDir 'SuppIndices/' parName '*']);
        parIndices = [];
        for thisRun = runIndices
            load([indicesDir 'SuppIndices/' parFiles(thisRun).name]);
            parIndices = [parIndices suppIndex];
        end
        analRow(1) = mean(parIndices); 
        
        %% get epoch durations for rivalry---------------------------------
        
        % get participant files
        parFiles = dir([indicesDir 'EpochDurations/' parName '*']);
        slowDurs = []; fastDurs = []; mixedDurs = [];
        for thisRun = runIndices
            load([indicesDir 'EpochDurations/' parFiles(thisRun).name]);
            slowDurs = [slowDurs durData(1,1)];
            fastDurs = [fastDurs durData(2,1)];
            mixedDurs = [mixedDurs durData(3,1)];
        end
        analRow(2) = mean(slowDurs);
        analRow(3) = mean(fastDurs);
        analRow(4) = mean(mixedDurs);
        
        %% get rivalry modulation indices----------------------------------
        
        % participants will not have rivalry modulation indices if SNR was 
        % too low (this was not an issue for the behavioral indices)
        if exist([indicesDir 'modIndices/' parName '.mat'], 'file')
            load([indicesDir 'modIndices/' parName '.mat']);
            analRow(5) = modIndex(1);
            analRow(6) = modIndex(2);
        end
        
         %% get transition rates-------------------------------------------
        
        % get participant files
        parFiles = dir([indicesDir 'tRate/' parName '*']);
        parIndices = [];
        for thisRun = runIndices
            load([indicesDir 'tRate/' parFiles(thisRun).name]);
            parIndices = [parIndices tRate];
        end
        analRow(7) = mean(parIndices); 
        
        analMatrix = [analMatrix; analRow];
    end
    save(['indices/marzRival_analMatrix']);
    R = corrcoef(analMatrix, 'rows', 'pairwise');
end

