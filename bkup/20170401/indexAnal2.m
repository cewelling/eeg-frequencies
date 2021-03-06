clearvars -except pValues kSize cutOff rValues

% load parameters
groupAnalysisParams

% Which type of run to analyze
analType = 'darts_dartRival'; %'darts'

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
% Column 22: Dom State length (Calculated from EEG transitions)
% Column 23: peak DomDur
% Column 24: low dom peak amp (peak to trough of itself)
% Column 25: high dom peak amp
% Column 26: low dom supp amp
% Column 27: high dom supp amp
% Column 28: low dom suppression index (supp amp / peak amp)
% Column 29: high dom suppression index
% Column 30: peak to trough period calculated from transition
% Column 31: peak to trough amplitude calculated from transition
% Column 32: oscillation rate calculated by taking fft of H-L
% Column 33: sim oscillation rate
% Column 34: normalized amplitude of fft at oscillation frequency

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
C{23} = 'peakDomDur';
C{24} = 'lowPkAmp';
C{25} = 'highPkAmp';
C{26} = 'lowSpAmp';
C{27} = 'highSpAmp';
C{28} = 'lSuppIndex';
C{29} = 'hSuppIndex';
C{30} = 'ptPeriod';
C{31} = 'ptAmp';
C{32} = 'ptRatio';
C{32} = 'OscFreq';
C{33} = 'OscFreqSIM';
C{34} = 'OscAmpNorm';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
analMatrix{1} = [];
analMatrix{2} = [];
ampMatrix = [];

for iGroup = 1:length(groupParNums)
    for parNum = groupParNums{iGroup}
        numFormat = '%02d';
        parName = [groupCodes{iGroup} num2str(parNum, numFormat)];
        analRow = nan(1, 34);
        
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
        
%         if altPeak
%             if exist(['peaks2/' parName '_' normType '.mat'], 'file')
%                 load(['peaks2/' parName '_' normType '.mat'])
%             else
%                 continue;
%             end
%         else
%             if exist(['peaks/' parName '_' normType '.mat'], 'file')
%                 load(['peaks/' parName '_' normType '.mat'])
%             else
%                 continue;
%             end
%         end
%         
%         % peak amps
%         dartRlowDom_amp = nanmax(nanmean(pLists{1}.f1, 1));
%         dartRhighDom_amp = nanmax(nanmean(pLists{2}.f2, 1));
%         dartSlowDom_amp = nanmax(nanmean(pLists{3}.f1, 1));
%         dartShighDom_amp = nanmax(nanmean(pLists{4}.f2, 1));
%         marzRlowDom_amp = nanmax(nanmean(pLists{5}.f1, 1));
%         marzRhighDom_amp = nanmax(nanmean(pLists{6}.f2, 1));
% %         
% %         ampRow = [dartRlowDom_amp dartRhighDom_amp dartSlowDom_amp dartShighDom_amp marzRlowDom_amp marzRhighDom_amp];
% %         ampMatrix = [ampMatrix; ampRow];
%         
%         if strcmp(analType, 'darts_dartRival')
%             % dominant amps (darts)
%             analRow(9) = dartRlowDom_amp;
%             analRow(10) = dartRhighDom_amp;
%             analRow(11) = mean([dartRlowDom_amp dartRhighDom_amp]);
%             
%             % suppression amps (darts)
%             lowDom_supp = nanmin(nanmean(pLists{1}.f2, 1));
%             highDom_supp = nanmin(nanmean(pLists{2}.f1, 1));
%         else
%             % dominant amps (marz)
%             analRow(9) = marzRlowDom_amp;
%             analRow(10) = marzRhighDom_amp;
%             analRow(11) = mean([marzRlowDom_amp marzRhighDom_amp]);
%             
%             % suppression amps (marz)
%             lowDom_supp = nanmin(nanmean(pLists{5}.f2, 1));
%             highDom_supp = nanmin(nanmean(pLists{6}.f1, 1));
%         end
%         
%         analRow(12) = lowDom_supp;
%         analRow(13) = highDom_supp;
%         analRow(14) = mean([lowDom_supp highDom_supp]);
%         
%         % domInd
%         if strcmp(analType, 'darts_dartRival')
%             domIndL = (dartRlowDom_amp - 1) / (dartSlowDom_amp - 1);
%             domIndH = (dartRhighDom_amp - 1) / (dartShighDom_amp - 1);
%             analRow(15) = domIndL;
%             analRow(16) = domIndH;
%             analRow(17) = mean([domIndL domIndH]);
%         end
%         
%         % dsInd
%         if strcmp(analType, 'darts_dartRival')
%             dsIndL = dartRlowDom_amp / lowDom_supp;
%             dsIndH = dartRhighDom_amp / highDom_supp;
%         else
%             dsIndL = marzRlowDom_amp / lowDom_supp;
%             dsIndH = marzRhighDom_amp / highDom_supp;
%         end
%         analRow(18) = dsIndL;
%         analRow(19) = dsIndH;
%         analRow(20) = mean([dsIndL dsIndH]);
        
        %% get periods from transition plots (time between high/low frequency overlaps closest to the transition overlap)
        
%         Pds = nan(1, 6); % 2 dartRival (h2l, l2h), 2 marzRival
%         
%         %load averaged transitions for the participant
%         %                 if exist(['transitions/' parName '_' normType '.mat'], 'file')
%         %                     load(['transitions/' parName '_' normType '.mat'])
%         %                 else
%         %                     continue;
%         %                 end
%         
%         if exist(['normTransitions/' parName '_' normType '.mat'], 'file')
%             load(['normTransitions/' parName '_' normType '.mat'])
%         else
%             continue;
%         end
%         tLists = ntLists;
%         
%         if strcmp(analType, 'darts')
%             tTypes = [1 2];
%         else
%             tTypes = [5 6]; % marz
%         end
%         
%         % find start and endpts of period, store pd length, dom state length
%         allPds = [];
%         domStates = [];
%         ptPds = [];
%         ptAmps = [];
%         ptRatio = [];
%         for tType = tTypes
%             if ~isempty(tLists{tType})
%                 
%                 %             %for averaging transitions before calculating period
%                 tLists{tType}.f1 = nanmean(tLists{tType}.f1, 1);
%                 tLists{tType}.f2 = nanmean(tLists{tType}.f2, 1);
%                 
%                 for t = 1:size(tLists{tType}.f1,1)
%                     
%                     f1trace = tLists{tType}.f1(t,:);
%                     f2trace = tLists{tType}.f2(t,:);
%                     
%                     % Normalize transitions
%                     
%                     % by lining up mean
%                     f1trace = f1trace - nanmean(f1trace);
%                     f2trace = f2trace - nanmean(f2trace);
%                     
%                     % by lining up min
%                     %                 f1trace = f1trace - min(f1trace);
%                     %                 f1trace = f1trace - min(f1trace);
%                     
%                     [leftCross, leftDist] = findCross(f1trace, f2trace, buttonPress, 'left');
%                     [rightCross, rightDist] = findCross(f1trace, f2trace, buttonPress, 'right');
%                     
%                     slopeInt = 50; % interval for determining pos vs neg slope
%                     
%                     if isnan(leftCross) || leftCross <= slopeInt
%                         transCross = rightCross;
%                         leftEnd = leftCross;
%                         rightEnd = findCross(f1trace, f2trace, transCross, 'right');
%                     elseif isnan(rightCross) || rightCross >= length(f1trace) - slopeInt
%                         transCross = leftCross;
%                         rightEnd = rightCross;
%                         [leftEnd, ~] = findCross(f1trace, f2trace, transCross, 'left');
%                     else                       
%                         % check in higher frequency (generally has larger dynamic range)
%                         f2leftDiff = f2trace(leftCross - slopeInt) - f2trace(leftCross + slopeInt);
%                         f2rightDiff = f2trace(rightCross - slopeInt) - f2trace(rightCross + slopeInt);
%                         
%                         if tType == 1
%                             if f2leftDiff <= f2rightDiff
%                                 transCross = leftCross;
%                                 rightEnd = rightCross;
%                                 [leftEnd, ~] = findCross(f1trace, f2trace, transCross, 'left');
%                             else
%                                 transCross = rightCross;
%                                 leftEnd = leftCross;
%                                 rightEnd = findCross(f1trace, f2trace, transCross, 'right');
%                             end
%                         elseif tType == 2
%                             if f2leftDiff >= f2rightDiff
%                                 transCross = leftCross;
%                                 rightEnd = rightCross;
%                                 [leftEnd, ~] = findCross(f1trace, f2trace, transCross, 'left');
%                             else
%                                 transCross = rightCross;
%                                 leftEnd = leftCross;
%                                 rightEnd = findCross(f1trace, f2trace, transCross, 'right');
%                             end
%                         end
%                     end
%                     
% %                     %identify the closer cross as the transition itself
% %                                         if rightDist < leftDist
% %                                             transCross = rightCross;
% %                                             leftEnd = leftCross;
% %                                             rightEnd = findCross(f1trace, f2trace, transCross, 'right');
% %                                             % if distances are equal, assume cross on the left is the
% %                                             % transition (there should be some lag)
% %                                         else
% %                     transCross = leftCross;
% %                     rightEnd = rightCross;
% %                     [leftEnd, ~] = findCross(f1trace, f2trace, transCross, 'left');
% %                                         end
%                     
%                     % store period
%                     Pds(tType) = rightEnd - leftEnd;
%                     allPds = [allPds Pds(tType)];
%                     
%                     % store domState length
%                     domStates = [domStates rightEnd - transCross transCross - leftEnd];
%                     
%                     % store peak-to-trough period
%                     
% %                     if isnan(leftEnd)|| isnan(rightEnd)
% %                         % initialize peak and trough indices
% %                         peakIs = nan(1,2);
% %                         troughIs = nan(1,2);
% %                         
% %                         % starting at crossing point, find the surrounding peaks and troughs
% %                         if tType == 1
% %                             upSteps = [-1 1]; % left (up), right (down) for f1
% %                         else
% %                             upSteps = [1 -1];
% %                         end
% %                         
% %                         traces = {f1trace f2trace};
% %                         
% %                         for freq = 1:2 %[{f1trace f1PeakI f1TroughI f1steps} {f2trace f2PeakI f2TroughI f2steps}]
% %                             % trace peak trough step
% %                             
% %                             % step up from transition cross
% %                             i = transCross;
% %                             step = upSteps(freq);
% %                             platCount = 0; % counts number points in plateau if top of peak is flat
% %                             found = 0;
% %                             while i > 1 && i < length(f1trace)
% %                                 [platCount, found] = stepUp(traces{freq}, i, platCount, step, found);
% %                                 if found == 1
% %                                     peakIs(freq) = i + round(platCount/2);
% %                                     break;
% %                                 end
% %                                 i = i+upSteps(freq);
% %                             end
% %                             
% %                             % step down from transition cross
% %                             i = transCross;
% %                             step = -upSteps(freq);
% %                             platCount = 0;
% %                             found = 0;
% %                             while i > 1 && i < length(f1trace)
% %                                 [platCount, found] = stepDown(traces{freq}, i, platCount, step, found);
% %                                 if found == 1
% %                                     troughIs(freq) = i - round(platCount/2);
% %                                     break;
% %                                 end
% %                                 i = i-upSteps(freq);
% %                             end
% %                         end
% %                         
% %                         f1PeakI = peakIs(1);
% %                         f1TroughI = troughIs(1);
% %                         f2PeakI = peakIs(2);
% %                         f2TroughI = troughIs(2);
% 
%                     if isnan(leftEnd)
%                         leftEnd = 1;
%                     end
%                     if isnan(rightEnd)
%                         rightEnd = length(f1trace);
%                     end
% %                     else
% %                         if isnan(transCross)
% %                             disp('hi')
% %                         end
%                         if tType == 1
%                             [~, f1PeakI] = nanmax(f1trace(leftEnd:transCross));
%                             f1PeakI = f1PeakI - 1 + leftEnd;
%                             [~, f1TroughI] = nanmin(f1trace(transCross:rightEnd));
%                             f1TroughI = f1TroughI - 1 + transCross;
%                             [~, f2PeakI] = nanmax(f2trace(transCross:rightEnd));
%                             f2PeakI = f2PeakI - 1 + transCross;
%                             [~, f2TroughI] = nanmin(f2trace(leftEnd:transCross));
%                             f2TroughI = f2TroughI - 1 + leftEnd;
%                             
%                         else
%                             [~, f1PeakI] = nanmax(f1trace(transCross:rightEnd));
%                             f1PeakI = f1PeakI - 1 + transCross;
%                             [~, f1TroughI] = nanmin(f1trace(leftEnd:transCross));
%                             f1TroughI = f1TroughI - 1 + leftEnd;
%                             [~, f2PeakI] = nanmax(f2trace(leftEnd:transCross));
%                             f2PeakI = f2PeakI - 1 + leftEnd;
%                             [~, f2TroughI] = nanmin(f2trace(transCross:rightEnd));
%                             f2TroughI = f2TroughI - 1 + transCross;
%                             
%                         end
% %                     end
% % if iGroup == 2
% %                     figure;
% %                     plot(f1trace)
% %                     hold on
% %                     plot(f2trace)
% %                     plot(f1PeakI, f1trace(f1PeakI), 'ro');
% %                     plot(f1TroughI, f1trace(f1TroughI), 'ro');
% %                     plot(f2PeakI, f2trace(f2PeakI), 'bo');
% %                     plot(f2TroughI, f2trace(f2TroughI), 'bo');
% %                     hold off
% % end
%                     
%                     ptPds = [ptPds abs(f1PeakI - f1TroughI) abs(f2PeakI - f2TroughI)];
%                     %ptPds = [ptPds abs(f2PeakI - f2TroughI)];
%                     %ptPds = [ptPds abs(peakIs(1) - troughIs(1)) abs(peakIs(2) - troughIs(2))];
%                     ptAmps = [ptAmps f2trace(f2PeakI) - f2trace(f2TroughI) f1trace(f1PeakI) - f1trace(f1TroughI)];
%                     ptRatio = [ptRatio f2trace(f2PeakI) - f2trace(f2TroughI) / f1trace(f1PeakI) - f1trace(f1TroughI)];
%                 end
%             end
%         end
%         
%         %if length(allPds) > minInstances
%         analRow(21) = nanmean(allPds);
%         %end
%         %if length(domStates) > minInstances
%         analRow(22) = nanmean(domStates);
%         %end
%         %if length(ptPds) > minInstances * 2
%         analRow(30) = nanmean(ptPds);
%         analRow(31) = nanmean(ptAmps);
%         analRow(32) = nanmean(ptRatio);
%         %end
        
        %% get domDur from peak plots
        
%         %load averaged transitions for the participant
%         %         if altPeak && exist(['peaks2/' parName '_' normType '.mat'], 'file')
%         %             load(['peaks2/' parName '_' normType '.mat'])
%         %         elseif ~altPeak && exist(['peaks/' parName '_' normType '.mat'], 'file')
%         %             load(['peaks/' parName '_' normType '.mat']);
%         %         else
%         %             continue;
%         %         end
%         
%         if exist(['normPeaks/' parName '_' normType '.mat'], 'file')
%             load(['normPeaks/' parName '_' normType '.mat'])
%         else
%             continue;
%         end
%         pLists = npLists;
%         
%         if strcmp(analType, 'darts')
%             pTypes = [1 2];
%         else
%             pTypes = [5 6]; % marz
%         end
%         
%         % find start and endpts of period, store pd length, dom state length
%         domStates = [];
%         lowPkAmp = []; % amplitude of aligned freq (peak to trough)
%         highPkAmp = [];
%         lowSpAmp = []; % amplitude of suppressed freq (peak to trough)
%         highSpAmp = [];
%         lSuppIndex = []; % amp of suppressed freq / amp of aligned freq
%         hSuppIndex = [];
%         
%         
%         for pType = pTypes
%             if ~isempty(pLists{pType})
%                 
%                 %             %TEMP for averaging transitions before calculating period
%                 %             tLists{tType}.f1 = nanmean(tLists{tType}.f1, 1);
%                 %             tLists{tType}.f2 = nanmean(tLists{tType}.f2, 1);
%                 
%                 for p = 1:size(pLists{pType}.f1,1)
%                     
%                     f1trace = pLists{pType}.f1(p,:);
%                     f2trace = pLists{pType}.f2(p,:);
%                     
%                     % Normalize peaks
%                     
%                     % by lining up mean
%                     f1trace = f1trace - nanmean(f1trace);
%                     f2trace = f2trace - nanmean(f2trace);
%                     %
%                     % by lining up min
%                     %                 f1trace = f1trace - min(f1trace);
%                     %                 f1trace = f1trace - min(f1trace);
%                     
%                     [leftEnd,~] = findCross(f1trace, f2trace, peakI, 'left');
%                     [rightEnd,~] = findCross(f1trace, f2trace, peakI, 'right');
%                     
%                     % store domState length
%                     domStates = [domStates rightEnd - leftEnd];
%                     
%                     % store peak amplitudes
%                     if pType == 1
%                         lowPkAmp = [lowPkAmp max(f1trace) - min(f1trace)];
%                         lowSpAmp = [lowSpAmp max(f2trace) - min(f2trace)];
%                         lSuppIndex = [lSuppIndex (lowSpAmp / lowPkAmp)];
%                     else
%                         highPkAmp = [highPkAmp max(f1trace) - min(f1trace)];
%                         highSpAmp = [highSpAmp max(f2trace) - min(f2trace)];
%                         hSuppIndex = [hSuppIndex (highSpAmp / highPkAmp)];
%                     end
%                 end
%             end
%         end
%         
%         if length(domStates) > minInstances % if one peak index has enough, all other peak indices should
%             analRow(23) = nanmean(domStates);
%             analRow(24) = nanmean(lowPkAmp);
%             analRow(25) = nanmean(highPkAmp);
%             analRow(26) = nanmean(lowSpAmp);
%             analRow(27) = nanmean(highSpAmp);
%             analRow(28) = nanmean(lSuppIndex);
%             analRow(29) = nanmean(hSuppIndex);
%         end
        
        %% get dart amps normalized to sim amps--------------------------------
        
        %         % No sim trials for marzipan, no modulation indices
        %         if strcmp(analType, 'marzipan')
        %             analRow(23) = NaN;
        %         else
        %             % participants will not have rivalry modulation indices if SNR was
        %             % too low (this was not an issue for the behavioral indices)
        %             if exist([indicesDir 'normAmps/' parName '.mat'], 'file')
        %                 load([indicesDir 'normAmps/' parName '.mat']);
        %                 analRow(23) = normAmp;
        %             else
        %                 analRow(23) = NaN;
        %             end
        %         end
        
        %% Load oscillation frequencies
        
        if exist([indicesDir 'rateFreqs/' parName '.mat'], 'file')
                load([indicesDir 'rateFreqs/' parName '.mat']);
                analRow(32) = oscFreq;
        end
            
        if exist([indicesDir 'rateFreqs/' parName '_SIM.mat'], 'file')
                load([indicesDir 'rateFreqs/' parName '_SIM.mat']);
                analRow(33) = oscFreq;
        end 
        
        if exist([indicesDir 'oscAmps/' parName '.mat'], 'file')
                load([indicesDir 'oscAmps/' parName '.mat']);
                rivAmp = oscAmp;
                
                if exist([indicesDir 'oscAmps/' parName '_SIM.mat'], 'file')
                    load([indicesDir 'oscAmps/' parName '_SIM.mat']);
                    normAmp = rivAmp / oscAmp;
                end
                analRow(34) = normAmp;
        end
        
        analMatrix{iGroup} = [analMatrix{iGroup}; analRow];
    end
    
    
    % calculate correlation coefficients
    [R, P] = corrcoef(analMatrix{iGroup}, 'rows', 'pairwise');
    %[R, P] = corr(analMatrix{iGroup}(:,5),analMatrix{iGroup}(:,32), 'type','spearman', 'rows', 'pairwise');

    
    
    % save full analysis matrices
    currMatrix = analMatrix{iGroup};
    if strcmp(analType, 'darts')
        save(['indices/' analGroupIDs{iGroup} 'darts_analMatrix2', 'currMatrix']);
    else
        save(['indices/' analGroupIDs{iGroup} 'marz_analMatrix2', 'currMatrix']);
    end
    
    % plot amps
%     figure
%     hold on
%     ax = {'DartRivalry-lowDom', 'DartRivalry-highDom', 'DartSim-lowDom', 'DartSim-highDom', 'MarzRivalry-lowDom', 'MarzRivalry-highDom'};
%     xaxis = {};
%     for i = 1:size(ampMatrix, 2)
%         xValues = i*ones(size(ampMatrix,1), 1);
%         scatter(xValues, ampMatrix(:,i));
%         xaxis = [xaxis ax{i}];
%     end
%     set(gca, 'XTick', [1:i], 'XTickLabel', xaxis)
    
    
    % plot potential correlations
    for i = [4:5] %[2 4:5] %size(fullAnalMatrix, 2)
        for j = [32] %i+1:size(fullAnalMatrix,2)
            
            %if abs(P(i,j)) < 0.1
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
                %text(xax(1) + (xax(2) - xax(1))/2, yax(1) + (yax(2) - yax(1))*.91,['y = ' num2str(p(1)) 'x + ' num2str(p(2))])
                line = p(1) .* analMatrix{iGroup}(validBoth,i) + p(2); % compute a new vector r that has matching datapoints in x
                plot(analMatrix{iGroup}(validBoth,i), line, '-'); %axis([-1 14 0 2])
                
                xlabel(['dartRival ' C{i}]);
                ylabel(['dartRival ' C{j}]);
                
                % save figures
                if strcmp(analType, 'darts_dartRival')
                    figName = [corrPlotDir label{iGroup} '_' C{i} 'vs' C{j} '_DARTS.jpg'];
                else
                    figName = [corrPlotDir label{iGroup} '_' C{i} 'vs' C{j} '_MARZ.jpg'];
                end
                saveas(gcf,figName,'jpg')
            %end
            
        end
     end
    %figure
    %hist(analMatrix{iGroup}(:,33))
    %title(analGroupIDs{iGroup})
end

% figure; 
% bar([nanmean(analMatrix{1}(:,30)) nanmean(analMatrix{2}(:,30))]);
% hold on
% ste1 = ste(analMatrix{1}(:,30));
% ste2 = ste(analMatrix{2}(:,30));
% errorbar([nanmean(analMatrix{1}(:,30)) nanmean(analMatrix{2}(:,30))], [ste1 ste2],'.');




