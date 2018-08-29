clearvars

% load parameters
groupAnalysisParams

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
fullAnalMatrix = [];
ampMatrix = [];

for iGroup = 1:length(groupParNums)
    for parNum = groupParNums{iGroup}
        numFormat = '%02d';
        parName = [groupCodes{iGroup} num2str(parNum, numFormat)];
        analRow = nan(1, 12);
          
        % get participant files
        parFiles = dir([bResultsDir parName '*']);
        
        %% get dartboard rivalry indices
        
        % dominance proportion     
        if exist([bResultsDir parName '_darts_propDom.txt'], 'file')
            analRow(1) = load([bResultsDir parName '_darts_propDom.txt']);
        end
        
        % avg dominant epoch duration
        if exist([bResultsDir parName '_darts_StateDurs.txt'])
            analRow(2) = load([bResultsDir parName '_darts_StateDurs.txt']);
        end
        
        % avg mixed epoch duration
        if exist([bResultsDir parName '_darts_UpDurs.txt'])
            analRow(3) = load([bResultsDir parName '_darts_UpDurs.txt']);
        end
        
        % Num transitions
        if exist([bResultsDir parName '_darts_Transitions.txt'])
            analRow(4) = load([bResultsDir parName '_darts_Transitions.txt']);
        end
        
        % Num switches
        if exist([bResultsDir parName '_darts_Switches.txt'])
            analRow(5) = load([bResultsDir parName '_darts_Switches.txt']);
        end
            
        %% get rivalry modulation indices----------------------------------
        
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
        
        %% get amplitudes (mean-->peak) and suppression amps (mean-->trough)
        if exist(['peaks/' parName '.mat'], 'file')
            load(['peaks/' parName '.mat']);
            
            % peak amps
            dartRlowDom_amp = nanmax(nanmean(pLists{1}.f1, 1));
            dartRhighDom_amp = nanmax(nanmean(pLists{2}.f2, 1));
            dartSlowDom_amp = nanmax(nanmean(pLists{3}.f1, 1));
            dartShighDom_amp = nanmax(nanmean(pLists{4}.f2, 1));
            marzRlowDom_amp = nanmax(nanmean(pLists{5}.f1, 1));
            marzShighDom_amp = nanmax(nanmean(pLists{6}.f2, 1));
            
            ampRow = [dartRlowDom_amp dartRhighDom_amp dartSlowDom_amp dartShighDom_amp marzRlowDom_amp marzShighDom_amp];
            ampMatrix = [ampMatrix; ampRow];
            
            analRow(8) = dartRlowDom_amp;
            analRow(9) = dartRhighDom_amp;
            
            % suppression amps
            dartRlowDom_supp = nanmin(nanmean(pLists{1}.f2, 1));
            dartRhighDom_supp = nanmin(nanmean(pLists{2}.f1, 1));
            
            analRow(10) = dartRlowDom_supp;
            analRow(11) = dartRhighDom_supp; 
        end
        
        %% get periods (time between high/low frequency overlaps closest to the transition overlap)
        
        Pds = nan(1, 6); % 2 dartRival (h2l, l2h), 2 marzRival
        
        % load averaged transitions for the participant
        if exist(['transitions/' parName '.mat'], 'file')
            load(['transitions/' parName '.mat'])
        end
        
        for tType = [1 2] %5 6] % rivalry cases only %1:length(tLists)
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
        analRow(12) = nanmean(Pds(1:2));
        %analRow(13) = nanmean(Pds(5:6));
        
        analMatrix{iGroup} = [analMatrix{iGroup}; analRow];
    end
    % create full analysis matrix
 fullAnalMatrix = [fullAnalMatrix; analMatrix{iGroup}]
end

% calculate correlation coefficients
[R, P] = corrcoef(fullAnalMatrix, 'rows', 'pairwise');

 save('indices/dartRival_analMatrix2');
 
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
for i = 1:size(fullAnalMatrix, 2)
    for j = i+1:size(fullAnalMatrix,2)
   
            if abs(P(i,j)) < 0.1
                figure
                hold on
                plot(analMatrix{1}(:,i), analMatrix{1}(:,j), 'ko')
                plot(analMatrix{2}(:,i), analMatrix{2}(:,j), 'bo')
                
                text(2,1,['r =' num2str(R(i,j)) '   p =' num2str(P(i,j))])               
                p = polyfit(fullAnalMatrix(:,i),fullAnalMatrix(:,j),1);   % p returns 2 coefficients fitting r = a_1 * x + a_2
                line = p(1) .* fullAnalMatrix(:,i) + p(2); % compute a new vector r that has matching datapoints in x     
                plot(fullAnalMatrix(:,i), line, '-'); %axis([-1 14 0 2])
                
                xlabel(C{i});
                ylabel(C{j});
            end
            
    end
end
 
    
    

