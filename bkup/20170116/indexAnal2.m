clearvars

% load parameters
groupAnalysisParams

% Which type of run to analyze
analType = 'darts'; %'darts'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the following matrix for index analysis:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Column 1: Dominance proportion
% Column 2: Dominant epoch durations
% Column 3: Mixed epoch durations
% Column 4: Number of transitions
% Column 5: Number of switches
% Column 6: Rivalry modulation index (low dom), peaks
% Column 7: Rivalry modulation index (high dom), peaks
% Column 8: Non-normalized amp of low dom peaks
% Column 9: Non-normalized amp of high dom peaks
% Column 10: Non-normalized suppression amp for low dom peaks
% Column 11: Non-normalized suppression amp for high dom peaks
% Column 12: dartRival period

C{1} = 'domProp';
C{2} = 'domDur';
C{3} = 'mixedDur';
C{4} = 'numTrans';
C{5} = 'numSwitches';
C{6} = 'lowDomMod';
C{7} = 'highDomMod';
C{8} = 'lowDomAmp';
C{9} = 'highDomAmp';
C{10} = 'lowDomSuppAmp';
C{11} = 'highDomSuppAmp';
C{12} = 'period';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
analMatrix{1} = [];
analMatrix{2} = [];
ampMatrix = [];

for iGroup = 1:length(groupParNums)
    for parNum = groupParNums{iGroup}
        numFormat = '%02d';
        parName = [groupCodes{iGroup} num2str(parNum, numFormat)];
        analRow = nan(1, 12);
        
        % get participant files
        parFiles = dir([bResultsDir parName '*']);
        
        %% get rivalry indices
        
        % dominance proportion
        if exist([bResultsDir parName '_' analType '_propDom.txt'], 'file')
            analRow(1) = load([bResultsDir parName '_' analType '_propDom.txt']);
        end
        
        % avg dominant epoch duration
        if exist([bResultsDir parName '_' analType '_StateDurs.txt'])
            analRow(2) = load([bResultsDir parName '_' analType '_StateDurs.txt']);
        end
        
        % avg mixed epoch duration
        if exist([bResultsDir parName '_' analType '_UpDurs.txt'])
            analRow(3) = load([bResultsDir parName '_' analType '_UpDurs.txt']);
        end
        
        % Num transitions
        if exist([bResultsDir parName '_' analType '_Transitions.txt'])
            analRow(4) = load([bResultsDir parName '_' analType '_Transitions.txt']);
        end
        
        % Num switches
        if exist([bResultsDir parName '_' analType '_Switches.txt'])
            analRow(5) = load([bResultsDir parName '_' analType '_Switches.txt']);
        end
        
        %% get rivalry modulation indices----------------------------------
        
        % No sim trials for marzipan, no modulation indices
        if strcmp(analType, 'marzipan')
            analRow(6) = NaN;
            analRow(7) = NaN;
        else
            % participants will not have rivalry modulation indices if SNR was
            % too low (this was not an issue for the behavioral indices)
            if exist([indicesDir 'modIndices/' parName '_lowDom.mat'], 'file')
                load([indicesDir 'modIndices/' parName '_lowDom.mat']);
                analRow(6) = lowDomMod;
            end
            
            if exist([indicesDir 'modIndices/' parName '_highDom.mat'], 'file')
                load([indicesDir 'modIndices/' parName '_highDom.mat']);
                analRow(7) = highDomMod;
            end
        end
        
        %% get amplitudes (mean-->peak) and suppression amps (mean-->trough)
        if exist(['peaks/' parName '.mat'], 'file')
            load(['peaks/' parName '.mat']);
            
            % peak amps
            dartRlowDom_amp = nanmax(nanmean(pLists{1}.f1, 1));
            dartRhighDom_amp = nanmax(nanmean(pLists{2}.f2, 1));
            dartSlowDom_amp = nanmax(nanmean(pLists{3}.f1, 1));
            dartShighDom_amp = nanmax(nanmean(pLists{4}.f2, 1));
            marzRlowDom_amp = nanmax(nanmean(pLists{5}.f1, 1));
            marzRhighDom_amp = nanmax(nanmean(pLists{6}.f2, 1));
            
            ampRow = [dartRlowDom_amp dartRhighDom_amp dartSlowDom_amp dartShighDom_amp marzRlowDom_amp marzRhighDom_amp];
            ampMatrix = [ampMatrix; ampRow];
            
            if strcmp(analType, 'darts')
                % dominant amps (darts)
                analRow(8) = dartRlowDom_amp;
                analRow(9) = dartRhighDom_amp;
                
                % suppression amps (darts)
                lowDom_supp = nanmin(nanmean(pLists{1}.f2, 1));
                highDom_supp = nanmin(nanmean(pLists{2}.f1, 1));
            else
                % dominant amps (marz)
                analRow(8) = marzRlowDom_amp;
                analRow(9) = marzRhighDom_amp;
                
                % suppression amps (marz)
                lowDom_supp = nanmin(nanmean(pLists{5}.f2, 1));
                highDom_supp = nanmin(nanmean(pLists{6}.f1, 1));
            end
            
            analRow(10) = lowDom_supp;
            analRow(11) = highDom_supp;
        end
        
        %% get periods (time between high/low frequency overlaps closest to the transition overlap)
        
        Pds = nan(1, 6); % 2 dartRival (h2l, l2h), 2 marzRival
        
        % load averaged transitions for the participant
        if exist(['transitions/' parName '.mat'], 'file')
            load(['transitions/' parName '.mat'])
        end
        
        if strcmp(analType, 'darts')
            tTypes = [1 2];
        else
            tTypes = [5 6]; % marz
        end
        
        for tType = tTypes
            [leftCross, leftDist] = findCross(tLists{tType}, buttonPress, 'left');
            [rightCross, rightDist] = findCross(tLists{tType}, buttonPress, 'right');
            
            % identify the closer cross as the transition itself
            if rightDist < leftDist
                transCross = rightCross;
                leftEnd = leftCross;
                rightEnd = findCross(tLists{tType}, transCross, 'right');
                % if distances are equal, assume cross on the left is the
                % transition (there should be some lag)
            else
                transCross = leftCross;
                rightEnd = rightCross;
                [leftEnd, ~] = findCross(tLists{tType}, transCross, 'left');
            end
            
            % store period
            Pds(tType) = rightEnd - leftEnd;
        end
        analRow(12) = nanmean(Pds(tType));
        
        analMatrix{iGroup} = [analMatrix{iGroup}; analRow];
    end
    
    % calculate correlation coefficients
    [R, P] = corrcoef(analMatrix{iGroup}, 'rows', 'pairwise');
    
    % save full analysis matrices
    currMatrix = analMatrix{iGroup};
    if strcmp(analType, 'darts')
        save(['indices/' analGroupIDs{iGroup} 'darts_analMatrix2', 'currMatrix']);
    else
        save(['indices/' analGroupIDs{iGroup} 'marz_analMatrix2', 'currMatrix']);
    end
    
    % plot amps
    figure
    hold on
    ax = {'DartRivalry-lowDom', 'DartRivalry-highDom', 'DartSim-lowDom', 'DartSim-highDom', 'MarzRivalry-lowDom', 'MarzRivalry-highDom'};
    xaxis = {};
    for i = 1:size(ampMatrix, 2)
        xValues = i*ones(size(ampMatrix,1), 1);
        scatter(xValues, ampMatrix(:,i));
        xaxis = [xaxis ax{i}];
    end
    set(gca, 'XTick', [1:i], 'XTickLabel', xaxis)
    
    % plot potential correlations
    for i = 1:5 %size(fullAnalMatrix, 2)
        for j = 6:12 %i+1:size(fullAnalMatrix,2)
            
            if abs(P(i,j)) < 0.1
                figure
                hold on
                plot(analMatrix{iGroup}(:,i), analMatrix{iGroup}(:,j), 'ko')
                title(analGroupIDs{iGroup})
                xax = xlim;
                yax = ylim;
                text(xax(1) + (xax(2) - xax(1))/2, yax(1) + (yax(2) - yax(1))*.95,['r =' num2str(R(i,j)) '   p =' num2str(P(i,j))])
                validi = ~isnan(analMatrix{iGroup}(:,i));
                validj = ~isnan(analMatrix{iGroup}(:,j));
                validBoth = validi & validj;
                p = polyfit(analMatrix{iGroup}(validBoth,i),analMatrix{iGroup}(validBoth,j),1);   % p returns 2 coefficients fitting r = a_1 * x + a_2
                line = p(1) .* analMatrix{iGroup}(validBoth,i) + p(2); % compute a new vector r that has matching datapoints in x
                plot(analMatrix{iGroup}(validBoth,i), line, '-'); %axis([-1 14 0 2])
                
                xlabel([analType ' ' C{i}]);
                ylabel([analType ' ' C{j}]);
                
                % save figures
                if strcmp(analType, 'darts')
                    figName = [corrPlotDir label{iGroup} '_' C{i} 'vs' C{j} '_DARTS.jpg'];
                else
                    figName = [corrPlotDir label{iGroup} '_' C{i} 'vs' C{j} '_MARZ.jpg'];
                end
                saveas(gcf,figName,'jpg')
            end
            
        end
    end
end
 
    
    

