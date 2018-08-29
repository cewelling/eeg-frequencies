function [lowSNR highSNR] = getIMsnrs( parName )

% load parameters
analysisParams

lowSNRs = [];
highSNRs = [];

for iRunType =  1 %1 % dartRival only
    cRunType = runTypes{iRunType};
    
    for runIndex = runIndices
        runName = [cRunType num2str(runIndex)];
        
        % load or (if not yet set) set common electrodes for each type of run based on SNR
        if strfind(cRunType, 'dart')
            if exist(['pre-processing/highSNRelecs/' parName '_darts_' num2str(commonElecs) 'of' num2str(numElecs) 'elecs.mat'], 'file')
                load(['pre-processing/highSNRelecs/' parName '_darts_' num2str(commonElecs) 'of' num2str(numElecs) 'elecs']);
            else
                elecs = setElecs(parName);
            end
        else
            if exist(['pre-processing/highSNRelecs/' parName '_marz_' num2str(commonElecs) 'of' num2str(numElecs) 'elecs.mat'], 'file')
                load(['pre-processing/highSNRelecs/' parName '_marz_' num2str(commonElecs) 'of' num2str(numElecs) 'elecs']);
            else
                elecs = setElecs(parName);
            end
        end
        for iTrial = 1:numTrials
            if exist(['pre-processing/highSNRelecs/' parName '_' runName '_trial' num2str(iTrial) '_' num2str(numElecs) 'elecs.mat'],'file')
                load(['pre-processing/highSNRelecs/' parName '_' runName '_trial' num2str(iTrial) '_' num2str(numElecs) 'elecs.mat']);
                if nanmean(maxSNRs(2,ismember(maxSNRs(1,:), elecs))) < minSNR
                    %if nanmean(lFreqSNRs(2,ismember(maxSNRs(1,:), elecs))) < 2 || ...
                    %        nanmean(hFreqSNRs(2,ismember(maxSNRs(1,:), elecs))) < 2
                    continue;
                else
                    indSNR = nanmean(maxSNRs(2,ismember(maxSNRs(1,:), elecs)));
                end
                load(['pre-processing/highSNRelecs/' parName '_' runName '_trial' num2str(iTrial) '_'  num2str(numElecs) 'elecs_IMs'])
                lowSNRs = [lowSNRs nanmean(lFreqSNRs(2,ismember(maxSNRs(1,:), elecs)))/indSNR];
                highSNRs = [highSNRs nanmean(hFreqSNRs(2, ismember(maxSNRs(1,:), elecs)))/indSNR];
            end
        end
    end
end

lowSNR = nanmean(lowSNRs);
highSNR = nanmean(highSNRs);

% saves the average of im SNR divided by average individual SNR for each trial
save([indicesDir 'imSNRs/' parName], 'lowSNR', 'highSNR'); 

end
