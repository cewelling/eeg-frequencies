% load group parameters
groupAnalysisParams

% load parameters
analysisParams

numFormat = '%02d'; % for participant number

for iGroup = 1:length(groupParNums)
    elecList = [];
    for parNum = groupParNums{iGroup};
        parName = [groupCodes{iGroup} num2str(parNum, numFormat)];
        
        % load electrodes used for analysis, add to electrode list
        % if file with electrodes does not exist, participant does not exist
        if exist(['highSNRelecs/' parName '_darts_' num2str(commonElecs) 'of' num2str(numElecs) 'elecs.mat'], 'file')
            load(['highSNRelecs/' parName '_darts_' num2str(commonElecs) 'of' num2str(numElecs) 'elecs']);            
            elecList = [elecList elecs];
        end
        
    end
    
    % make histogram
    figure
    hist(elecList)
    title(analGroupIDs{iGroup});
end