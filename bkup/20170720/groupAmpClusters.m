% For an entire group, do amplitudes cluster according to perceptual report? 
% (Dominant vs. Mixed vs. Suppressed?)

% load parameters
groupAnalysisParams

for iRunType = 1:length(runTypes) % just input 1 for dartRival only
    cRunType = runTypes{iRunType};
    
    for iGroup = 1:length(groupParNums)
        
        % space to store group amplitudes corresponding to each percept
        hMixedG = [];
        hCorrG = []; % corresponding percept (high frequency if its the high frequency trace)
        hOppG = []; % opposite percept (low frequency of its the high frequency trace)
        lMixedG = [];
        lCorrG = [];
        lOppG = [];
        
        for parNum = groupParNums{iGroup}
            numFormat = '%02d';
            parName = [groupCodes{iGroup} num2str(parNum, numFormat)];
            
            % load amplitudes
            if exist(['ampClusters/' parName '_' cRunType '.mat'], 'file')
                load(['ampClusters/' parName '_' cRunType '.mat'])
            else
                continue;
            end
            
            % collect amplitudes
            hMixedG = [hMixedG nanmean(hMixedAmps)];
            hCorrG = [hCorrG nanmean(hCorrAmps)];
            hOppG = [hOppG nanmean(hOppAmps)];
            lMixedG = [lMixedG nanmean(lMixedAmps)];
            lCorrG = [lCorrG nanmean(lCorrAmps)];
            lOppG = [lOppG nanmean(lOppAmps)];
        end
        
        % Create group plots
        
        % Plot amplitude clusters for low frequency band
        figure
        bar([1 2 3], [nanmean(lCorrG) nanmean(lMixedG) nanmean(lOppG)]);
        hold on
        scatterX = [ones(1, length(lCorrG)) (ones(1, length(lMixedG)) + 1) (ones(1, length(lOppG)) + 2)];
        scatter(scatterX, [lCorrG lMixedG lOppG]);
        errorbar([1 2 3], [nanmean(lCorrG) nanmean(lMixedG) nanmean(lOppG)], [ste(lCorrG) ste(lMixedG) ste(lOppG)], 'k.','LineWidth',2);
        title([analGroupIDs{iGroup} ' ' cRunType ' : low frequency band'])
        xlabel('     Dominant                                       Mixed                                     Suppressed   ')
        
        % Plot amplitude clusters for high frequency band
        figure
        bar([1 2 3], [nanmean(hCorrG) nanmean(hMixedG) nanmean(hOppG)]);
        hold on
        scatterX = [ones(1, length(hCorrG)) (ones(1, length(hMixedG)) + 1) (ones(1, length(hOppG)) + 2)];
        scatter(scatterX, [hCorrG hMixedG hOppG]);
        errorbar([1 2 3], [nanmean(hCorrG) nanmean(hMixedG) nanmean(hOppG)], [ste(hCorrG) ste(hMixedG) ste(hOppG)], 'k.','LineWidth',2);
        title([analGroupIDs{iGroup} ' ' cRunType ' : high frequency band'])
        xlabel('     Dominant                                       Mixed                                     Suppressed   ')
        
        %         histEdges = -3:.2:3;
        %         figure
        %         subplot(3,1,1)
        %         histogram(lCorrG, histEdges)
        %         title([analGroupIDs{iGroup} ': low frequency band'])
        %         subplot(3,1,2)
        %         histogram(lMixedG, histEdges)
        %         subplot(3,1,3)
        %         histogram(lOppG, histEdges)
        %
        %         figure
        %         subplot(3,1,1)
        %         histogram(hCorrG, histEdges)
        %         title([analGroupIDs{iGroup} ': high frequency band'])
        %         subplot(3,1,2)
        %         histogram(hMixedG, histEdges)
        %         subplot(3,1,3)
        %         histogram(hOppG, histEdges)
    end
end