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
% Column 8: Avg of 6 & 7
% Column 9: Non-normalized amp of low dom peaks
% Column 10: Non-normalized amp of high dom peaks
% Column 11: Avg of 8 & 9
% Column 12: Non-normalized suppression amp for low dom peaks
% Column 13: Non-normalized suppression amp for high dom peaks
% Column 14: Avg of 11 & 12
% Column 15: low dom dominance index (riv amp / sim amp)
% Column 16: high dom dominance index (riv amp / sim amp)
% Column 17: Avg of 14 & 15
% Column 18: low dom ds index (riv amp / supp amp)
% Column 19: high dom ds index (riv amp / supp amp)
% Column 20: Avg of 17 % 18
% Column 21: dartRival period
% Column 22: Dom State length (Calculated from EEG data)
% Column 23: normAmps (dart amps normalized to sim amps)

C{1} = 'domProp';
C{2} = 'domDur';
C{3} = 'mixedDur';
C{4} = 'numTrans';
C{5} = 'numSwitches';
C{6} = 'lowDomMod';
C{7} = 'highDomMod';
C{8} = 'meanMod';
C{9} = 'lowDomAmp';
C{10} = 'highDomAmp';
C{11} = 'meanAmp';
C{12} = 'lowDomSuppAmp';
C{13} = 'highDomSuppAmp';
C{14} = 'meanSuppAmp';
C{15} = 'lowDomInd';
C{16} = 'highDomInd';
C{17} = 'meanDomInd';
C{18} = 'lowDSind';
C{19} = 'highDSind';
C{20} = 'meanDSind';
C{21} = 'period';
C{22} = 'domDurEEG';
C{23} = 'normAmp';

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
        else
            continue;
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
                analRow(8) = mean([lowDomMod highDomMod]);
            end          
        end
        
        %% get amplitudes (mean-->peak) and suppression amps (mean-->trough), calculate domInd and dsInd
        
        if altPeak
            if strcmp(normType, 'none')
                if exist(['peaksNoNorm2/' parName '.mat'], 'file')
                    load(['peaksNoNorm2/' parName '.mat'])
                else
                    continue;
                end
            elseif exist(['peaks2/' parName '.mat'], 'file')
                load(['peaks2/' parName '.mat'])
            else
                continue;
            end
        elseif strcmp(normType, 'none')
            if exist(['peaksNoNorm/' parName '.mat'], 'file')
                load(['peaksNoNorm/' parName '.mat'])
            else
                continue;
            end
        else
            if exist(['peaks/' parName '.mat'], 'file')
                load(['peaks/' parName '.mat'])
            else
                continue;
            end
        end
        
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
            analRow(9) = dartRlowDom_amp;
            analRow(10) = dartRhighDom_amp;
            analRow(11) = mean([dartRlowDom_amp dartRhighDom_amp]);
            
            % suppression amps (darts)
            lowDom_supp = nanmin(nanmean(pLists{1}.f2, 1));
            highDom_supp = nanmin(nanmean(pLists{2}.f1, 1));
        else
            % dominant amps (marz)
            analRow(9) = marzRlowDom_amp;
            analRow(10) = marzRhighDom_amp;
            analRow(11) = mean([marzRlowDom_amp marzRhighDom_amp]);
            
            % suppression amps (marz)
            lowDom_supp = nanmin(nanmean(pLists{5}.f2, 1));
            highDom_supp = nanmin(nanmean(pLists{6}.f1, 1));
        end
        
        analRow(12) = lowDom_supp;
        analRow(13) = highDom_supp;
        analRow(14) = mean([lowDom_supp highDom_supp]);
        
        % domInd
        if strcmp(analType, 'darts')
            domIndL = (dartRlowDom_amp - 1) / (dartSlowDom_amp - 1);
            domIndH = (dartRhighDom_amp - 1) / (dartShighDom_amp - 1);
            analRow(15) = domIndL;
            analRow(16) = domIndH;
            analRow(17) = mean([domIndL domIndH]);
        end
        
        % dsInd
        if strcmp(analType, 'darts')
            dsIndL = dartRlowDom_amp / lowDom_supp;
            dsIndH = dartRhighDom_amp / highDom_supp;
        else
            dsIndL = marzRlowDom_amp / lowDom_supp;
            dsIndH = marzRhighDom_amp / highDom_supp;
        end
        analRow(18) = dsIndL;
        analRow(19) = dsIndH;
        analRow(20) = mean([dsIndL dsIndH]);
        
        %% get periods (time between high/low frequency overlaps closest to the transition overlap)
        
        Pds = nan(1, 6); % 2 dartRival (h2l, l2h), 2 marzRival
        
        % load averaged transitions for the participant
        if allTrans
            if strcmp(normType, 'none')
                if exist(['transNoNorm2/' parName '.mat'], 'file')
                    load(['transNoNorm2/' parName '.mat'])
                else
                    continue;
                end
            elseif exist(['transitions2/' parName '.mat'], 'file')
                load(['transitions2/' parName '.mat'])
            else
                continue;
            end
        elseif strcmp(normType, 'none')
            if exist(['transNoNorm/' parName '.mat'], 'file')
                load(['transNoNorm/' parName '.mat'])
            else
                continue;
            end
        else
            if exist(['transitions/' parName '.mat'], 'file')
                load(['transitions/' parName '.mat'])
            else
                continue;
            end
        end
        
        if strcmp(analType, 'darts')
            tTypes = [1 2];
        else
            tTypes = [5 6]; % marz
        end
        
        % find start and endpts of period, store pd length, dom state length
        allPds = [];
        domStates = [];
        for tType = tTypes
            for t = 1:size(tLists{tType}.f1,1)
                
                f1trace = tLists{tType}.f1(t,:);
                f2trace = tLists{tType}.f2(t,:);
                
                % Normalize transitions
                
                % by lining up mean
                f1trace = f1trace - mean(f1trace);
                f2trace = f2trace - mean(f2trace);
                
                % by lining up min
%                 f1trace = f1trace - min(f1trace);
%                 f1trace = f1trace - min(f1trace);
                
                [leftCross, leftDist] = findCross(f1trace, f2trace, buttonPress, 'left');
                [rightCross, rightDist] = findCross(f1trace, f2trace, buttonPress, 'right');
                
                % identify the closer cross as the transition itself
                if rightDist < leftDist
                    transCross = rightCross;
                    leftEnd = leftCross;
                    rightEnd = findCross(f1trace, f2trace, transCross, 'right');
                    % if distances are equal, assume cross on the left is the
                    % transition (there should be some lag)
                else
                    transCross = leftCross;
                    rightEnd = rightCross;
                    [leftEnd, ~] = findCross(f1trace, f2trace, transCross, 'left');
                end
                
                % store period
                Pds(tType) = rightEnd - leftEnd;
                allPds = [allPds Pds(tType)];
                
                % store domState length
                domStates = [domStates rightEnd - transCross transCross - leftEnd];
            end
        end
        analRow(21) = nanmean(allPds);
        analRow(22) = nanmean(domStates);
        
        %% get dart amps normalized to sim amps--------------------------------
        
        % No sim trials for marzipan, no modulation indices
        if strcmp(analType, 'marzipan')
            analRow(23) = NaN;
        else
            % participants will not have rivalry modulation indices if SNR was
            % too low (this was not an issue for the behavioral indices)
            if exist([indicesDir 'normAmps/' parName '.mat'], 'file')
                load([indicesDir 'normAmps/' parName '.mat']);
                analRow(23) = normAmp;
            else
                analRow(23) = NaN;
            end
        end
        
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
    for i = 1:3 %size(fullAnalMatrix, 2)
        for j = [21 22 23] %i+1:size(fullAnalMatrix,2)
            
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
 
    
    

