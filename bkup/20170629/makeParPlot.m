function [ tLists, tLabels, pLists, pLabels] = makeParPlot( parName, iPar, date)
% makeParPlot(parName, date) plots averaged transitions and peaks for each
% participant. Transitions and peaks with SNR < minSNR and instances <
% minInstances are excluded in this step. Lists of transitions and peaks
% are saved in the "transitions" and "peaks" folders, respectively.
% IMPORTANT: If planning to make group plots for more than one type of run
% (e.g. rivalry and simulation), makeParPlot must be run with both
% dartRival AND dartSIM chosen AT THE SAME TIME in analysisParams.m.
%
% [tLists, tLabels, pLists, pLabels] = makeParPlot(parName, date) returns
% lists of transitions, peaks, and their corresponding labels. tLists is a
% cell array, with each cell containing a struct corresponding to a type of
% transition as follows:
% 1. dartRival low to high
% 2. dartRival high to low
% 3. dartRival low to low
% 4. dartRival high to high
% 5. dartSim low to high
% 6. dartSim high to low
% 7. marzRival low to high
% 8. marzRival high to low
%
% Similarly, pLists is a cell array, with each cell containing a struct
% corresponding to a type of peak as follows:
% 1. dartRival low frequency peak
% 2. dartRival high frequency peak
% 3. dartSim low frequency peak
% 4. dartSim high frequency peak
% 5. marzRival low frequency peak
% 6. marzRival high frequency peak
%
% Each struct then has two fields, f1 and f2, that contain lists of the low
% and high frequency traces for that type of transition or peak.
%
% Called from: analysisController.m
% Dependencies: analysisParams.m, aggTransitions.m, aggPeaks.m

%% set-up

%import parameters
analysisParams

% transition case codes and peak types
l2h = 1; %low to high
h2l = 2;
l2l = 3;
h2h = 4;

slow = 1; % peaks in slow freq series
fast = 2;

analElecs = electrodeSet.nums;

for iElec = [1:length(analElecs) 0] % 0 for the average of the selected electrodes
    if iElec ~= 0 % for individual channels
        thisElec = analElecs(iElec);
    end
    
    % allocate space to store transitions
    l2hRival.f1 = []; l2hRival.f2 = []; l2hRival.tPoints = []; l2hRival.tTrials = [];
    h2lRival.f1 = []; h2lRival.f2 = []; h2lRival.tPoints = []; h2lRival.tTrials = [];
    l2lRival.f1 = []; l2lRival.f2 = []; l2lRival.tPoints = []; l2lRival.tTrials = [];
    h2hRival.f1 = []; h2hRival.f2 = []; h2hRival.tPoints = []; h2hRival.tTrials = [];
    l2hSim.f1 = []; l2hSim.f2 = []; l2hSim.tPoints = []; l2hSim.tTrials = [];
    h2lSim.f1 = []; h2lSim.f2 = []; h2lSim.tPoints = []; h2lSim.tTrials = [];
    l2hMarzRival.f1 = []; l2hMarzRival.f2 = []; l2hMarzRival.tPoints = []; l2hMarzRival.tTrials = [];
    h2lMarzRival.f1 = []; h2lMarzRival.f2 = []; h2lMarzRival.tPoints = []; h2lMarzRival.tTrials = [];
    l2lMarzRival.f1 = []; l2lMarzRival.f2 = []; l2lMarzRival.tPoints = []; l2lMarzRival.tTrials = [];
    h2hMarzRival.f1 = []; h2hMarzRival.f2 = []; h2hMarzRival.tPoints = []; h2hMarzRival.tTrials = [];
    
    % allocate space to store peaks
    slowPkRival.f1 = []; slowPkRival.f2 = [];
    fastPkRival.f1 = []; fastPkRival.f2 = [];
    slowPkSim.f1 = []; slowPkSim.f2 = [];
    fastPkSim.f1 = []; fastPkSim.f2 = [];
    slowPkMarzRival.f1 = []; slowPkMarzRival.f2 = [];
    fastPkMarzRival.f1 = []; fastPkMarzRival.f2 = [];
    
    %% collect all transitions and peaks for each participant
    
    for iRunType = 1:length(runTypes)
        cRunType = runTypes{iRunType};
        
        runStart = 1; % keep track of where current run starts in storage matrices
        
        for runIndex = runIndices
            runName = [cRunType num2str(runIndex)];
            
            % Get the appropriate EEG file for this run
            EEGfile = [eegDir parName '_' runName '_' date '.bdf'];
            
            % skip over participants / runs that don't exist
            if ~exist(EEGfile, 'file')
                continue;
            end
            
            % prepare to load transitions for this run
            if strcmp(cRunType, 'dartSim') && useSimSchedule
                tFolder = 'transitions/tListsRLS_simSched/';
            else
                tFolder = 'transitions/tListsRLS/';
            end
            
            if iElec ~= 0 % individual electrodes
                % load transitions for this run
                if exist([tFolder parName '_' runName '_' normType '_' analFreqLabel '_' num2str(thisElec) '.mat'], 'file') 
                    load([tFolder parName '_' runName '_' normType '_' analFreqLabel '_' num2str(thisElec) '.mat'])
                else
                    aggTransitions(parName, iPar, runName, date, EEGfile);
                    load([tFolder parName '_' runName '_' normType '_' analFreqLabel '_' num2str(thisElec) '.mat'])
                end           
            else % electrodes averaged together
                % load transitions for this run
                if exist([tFolder parName '_' runName '_' normType '_' analFreqLabel '_' electrodeSet.name '.mat'], 'file')
                    load([tFolder parName '_' runName '_' normType '_' analFreqLabel '_' electrodeSet.name '.mat'])
                else
                    transitions = aggTransitions(parName, iPar, runName, date, EEGfile);
                end
                
                % load peaks for this run
                if altPeak &&  ~exist(['peaks/altpListsRLS/' parName '_' runName '_' normType '.mat'], 'file') || ...
                        ~altPeak && ~exist(['peaks/pListsRLS/' parName '_' runName '_' normType '.mat'], 'file')
                    peaks = aggPeaks(parName, runName, date, EEGfile);
                elseif altPeak &&  exist(['peaks/altpListsRLS/' parName '_' runName '_' normType '.mat'], 'file')
                    load(['peaks/altpListsRLS/' parName '_' runName '_' normType '.mat'])
                else
                    load(['peaks/pListsRLS/' parName '_' runName '_' normType '.mat'])
                end
            end
            
            % determine which transitions fit criteria
            tKeep = [];
            for cCode = 1:4 % l2h, h2l, l2l, h2h
                criteria = transitions(cCode).crit;
                tKeep{cCode} = [];
                if isempty(criteria)
                    continue;
                end
                for i = 1:size(criteria, 1)
                    if criteria(i,7) >= minSNR && criteria(i,8) >= minSNR % exclude transitions by SNR (both frequencies must be above threshold)
                        if criteria(i,1) >= dom1Min && criteria(i,1) <= dom1Max && criteria(i,4) <= gapMax
                            if criteria(i,2) >= dom2Min && criteria(i,2) < dom2Max
                                if ~noMixed && ((criteria(i,3) >= mixedMin && criteria(i,3) <= mixedMax))
                                    tKeep{cCode} = [tKeep{cCode} i];
                                elseif noMixed && criteria(i,3) == 0
                                    tKeep{cCode} = [tKeep{cCode} i];
                                elseif iRunType == 2 % simulation trials, no mixed states
                                    tKeep{cCode} = [tKeep{cCode} i];
                                end
                            end
                        end
                    end
                end
            end
            
            if ~iElec
                % determine which peaks fit criteria
                pKeep = [];
                for domFreq = 1:2 % low, high
                    peakDurs = peaks(domFreq).durs;
                    peakSNRs = peaks(domFreq).snrs;
                    
                    pKeep{domFreq} = [];
                    if isempty(peakDurs)
                        continue;
                    end
                    for i = 1:size(peakDurs, 1)
                        if peakDurs(i) > minPwidth && peakDurs(i) < maxPwidth % exclude peaks by SNR
                            if peakSNRs(i) >= minSNR
                                pKeep{domFreq} = [pKeep{domFreq} i];
                            end
                        end
                    end
                end
            end
            
            if strfind(runName, 'dartRival')
                % dartboard rival transitions
                l2hRival.f1 = [l2hRival.f1; transitions(l2h).f1(tKeep{1}, :)];
                l2hRival.f2 = [l2hRival.f2; transitions(l2h).f2(tKeep{1}, :)];
                if ~isempty(transitions(l2h).crit)
                    l2hRival.tPoints = [l2hRival.tPoints; transitions(l2h).crit(tKeep{1}, 9)];
                    l2hRival.tTrials = [l2hRival.tTrials; transitions(l2h).crit(tKeep{1}, 10)];
                end
                h2lRival.f1 = [h2lRival.f1; transitions(h2l).f1(tKeep{2}, :)];
                h2lRival.f2 = [h2lRival.f2; transitions(h2l).f2(tKeep{2}, :)];
                if ~isempty(transitions(h2l).crit)
                    h2lRival.tPoints = [h2lRival.tPoints; transitions(h2l).crit(tKeep{2}, 9)];
                    h2lRival.tTrials = [h2lRival.tTrials; transitions(h2l).crit(tKeep{2}, 10)];
                end
                l2lRival.f1 = [l2lRival.f1; transitions(l2l).f1(tKeep{3}, :)];
                l2lRival.f2 = [l2lRival.f2; transitions(l2l).f2(tKeep{3}, :)];
                if ~isempty(transitions(l2l).crit)
                    l2lRival.tPoints = [l2lRival.tPoints; transitions(l2l).crit(tKeep{3}, 9)];
                    l2lRival.tTrials = [l2lRival.tTrials; transitions(l2l).crit(tKeep{3}, 10)];
                end
                h2hRival.f1 = [h2hRival.f1; transitions(h2h).f1(tKeep{4}, :)];
                h2hRival.f2 = [h2hRival.f2; transitions(h2h).f2(tKeep{4}, :)];
                if ~isempty(transitions(h2h).crit)
                    h2hRival.tPoints = [h2hRival.tPoints; transitions(h2h).crit(tKeep{4}, 9)];
                    h2hRival.tTrials = [h2hRival.tTrials; transitions(h2h).crit(tKeep{4}, 10)];
                end
                
                % dartboard rival peaks
                if ~iElec
                    slowPkRival.f1 = [slowPkRival.f1; peaks(slow).f1(pKeep{1}, :)];
                    slowPkRival.f2 = [slowPkRival.f2; peaks(slow).f2(pKeep{1}, :)];
                    fastPkRival.f1 = [fastPkRival.f1; peaks(fast).f1(pKeep{2}, :)];
                    fastPkRival.f2 = [fastPkRival.f2; peaks(fast).f2(pKeep{2}, :)];
                end
                
            elseif strfind(runName, 'dartSim')
                % dartboard sim transitions
                l2hSim.f1 = [l2hSim.f1; transitions(l2h).f1(tKeep{1}, :)];
                l2hSim.f2 = [l2hSim.f2; transitions(l2h).f2(tKeep{1}, :)];
                if ~isempty(transitions(l2h).crit)
                    l2hSim.tPoints = [l2hSim.tPoints; transitions(l2h).crit(tKeep{1}, 9)];
                    l2hSim.tTrials = [l2hSim.tTrials; transitions(l2h).crit(tKeep{1}, 10)];
                end
                h2lSim.f1 = [h2lSim.f1; transitions(h2l).f1(tKeep{2}, :)];
                h2lSim.f2 = [h2lSim.f2; transitions(h2l).f2(tKeep{2}, :)];
                if ~isempty(transitions(h2l).crit)
                    h2lSim.tPoints = [h2lSim.tPoints; transitions(h2l).crit(tKeep{2}, 9)];
                    h2lSim.tTrials = [h2lSim.tTrials; transitions(h2l).crit(tKeep{2}, 10)];
                end
                
                % dartboard sim peaks
                if ~iElec
                    slowPkSim.f1 = [slowPkSim.f1; peaks(slow).f1(pKeep{1}, :)];
                    slowPkSim.f2 = [slowPkSim.f2; peaks(slow).f2(pKeep{1}, :)];
                    fastPkSim.f1 = [fastPkSim.f1; peaks(fast).f1(pKeep{2}, :)];
                    fastPkSim.f2 = [fastPkSim.f2; peaks(fast).f2(pKeep{2}, :)];
                end
            
            elseif strfind(runName, 'marzRival')
                if ~iElec
                    % marz rival transitions
                    l2hMarzRival.f1 = [l2hMarzRival.f1; transitions(l2h).f1(tKeep{1}, :)];
                    l2hMarzRival.f2 = [l2hMarzRival.f2; transitions(l2h).f2(tKeep{1}, :)];
                    if ~isempty(transitions(l2h).crit)
                        l2hMarzRival.tPoints = [l2hMarzRival.tPoints; transitions(l2h).crit(tKeep{1}, 9)];
                        l2hMarzRival.tTrials = [l2hMarzRival.tTrials; transitions(l2h).crit(tKeep{1}, 10)];
                    end
                    h2lMarzRival.f1 = [h2lMarzRival.f1; transitions(h2l).f1(tKeep{2}, :)];
                    h2lMarzRival.f2 = [h2lMarzRival.f2; transitions(h2l).f2(tKeep{2}, :)];
                    if ~isempty(transitions(h2l).crit)
                        h2lMarzRival.tPoints = [h2lMarzRival.tPoints; transitions(h2l).crit(tKeep{1}, 9)];
                        h2lMarzRival.tTrials = [h2lMarzRival.tTrials; transitions(h2l).crit(tKeep{1}, 10)];
                    end
                    l2lMarzRival.f1 = [l2lMarzRival.f1; transitions(l2l).f1(tKeep{3}, :)];
                    l2lMarzRival.f2 = [l2lMarzRival.f2; transitions(l2l).f2(tKeep{3}, :)];
                    if ~isempty(transitions(l2l).crit)
                        l2lMarzRival.tPoints = [l2lMarzRival.tPoints; transitions(l2l).crit(tKeep{1}, 9)];
                        l2lMarzRival.tTrials = [l2lMarzRival.tTrials; transitions(l2l).crit(tKeep{1}, 10)];
                    end
                    h2hMarzRival.f1 = [h2hMarzRival.f1; transitions(h2h).f1(tKeep{4}, :)];
                    h2hMarzRival.f2 = [h2hMarzRival.f2; transitions(h2h).f2(tKeep{4}, :)];
                    if ~isempty(transitions(h2h).crit)
                        h2hMarzRival.tPoints = [h2hMarzRival.tPoints; transitions(h2h).crit(tKeep{1}, 9)];
                        h2hMarzRival.tTrials = [h2hMarzRival.tTrials; transitions(h2h).crit(tKeep{1}, 10)];
                    end
                    
                    % marz rival peaks
                    slowPkMarzRival.f1 = [slowPkMarzRival.f1; peaks(slow).f1(pKeep{1}, :)];
                    slowPkMarzRival.f2 = [slowPkMarzRival.f2; peaks(slow).f2(pKeep{1}, :)];
                    fastPkMarzRival.f1 = [fastPkMarzRival.f1; peaks(fast).f1(pKeep{2}, :)];
                    fastPkMarzRival.f2 = [fastPkMarzRival.f2; peaks(fast).f2(pKeep{2}, :)];
                end
            end
            
            if iElec ~= 0 % save by run for individual channels
                % take subset of transitions from just this run, prepare to save
                l2hRival_run = structSubset(l2hRival, runStart);
                h2lRival_run = structSubset(h2lRival, runStart);
                l2lRival_run = structSubset(l2lRival, runStart);
                h2hRival_run = structSubset(h2hRival, runStart);
                l2hSim_run = structSubset(l2hSim, runStart);
                h2lSim_run = structSubset(h2lSim, runStart);
%                 l2hMarzRival_run = structSubset(l2hMarzRival, runStart);
%                 h2lMarzRival_run = structSubset(h2lMarzRival, runStart);
%                 l2lMarzRival_run = structSubset(l2lMarzRival, runStart);
%                 h2hMarzRival_run = structSubset(h2hMarzRival, runStart);
                
                % save transitions by run, just l2h and h2l (for svm)
                if ~isempty(strfind(runName, 'dartRival'))
                    run_tLists = {l2hRival_run h2lRival_run};
                    runStart = runStart + size(l2hRival_run.f1, 1); % update with start index of next run (in transition storage matrices)
                elseif ~isempty(strfind(runName, 'dartSim'))
                    run_tLists = {l2hSim_run h2lSim_run};
                    runStart = runStart + size(l2hSim_run.f1, 1); % update with start index of next run (in transition storage matrices)
                end
                %run_tLists = {l2hRival_run h2lRival_run l2lRival_run h2hRival_run l2hSim_run h2lSim_run l2hMarzRival_run h2lMarzRival_run l2lMarzRival_run h2hMarzRival_run};
                
                save(['transitions/transitionsByRun/' parName '_' runName '_' normType '_' analFreqLabel '_' num2str(thisElec)], 'run_tLists');          
            end
        end
    end
end

tLists = {l2hRival h2lRival l2lRival h2hRival l2hSim h2lSim l2hMarzRival h2lMarzRival l2lMarzRival h2hMarzRival};

%% Normalize rivalry transitions to simulation transitions

% Find average sim amplitude
f1simAmps = [];
f2simAmps = [];
f1PTerror = [];
f2PTerror = [];
for tType = 5:6 % 5: sim low to high; 6: sim high to low
    
    % average transitions
    f1trace = nanmean(tLists{tType}.f1,1);
    f2trace = nanmean(tLists{tType}.f2,1);
    
    % error across transitions
    f1error = ste(tLists{tType}.f1);
    f2error = ste(tLists{tType}.f2);
    
    if ~isnan(f1trace) % no f1 transitions means no f2 transitions either
        
        % find peak and trough of averaged transition
        
        % demean time courses
        f1trace = f1trace - nanmean(f1trace);
        f2trace = f2trace - nanmean(f2trace);
        
        % identify where they cross on both sides of button press
        [leftCross, ~] = findCross(f1trace, f2trace, buttonPress, 'left');
        [rightCross, ~] = findCross(f1trace, f2trace, buttonPress, 'right');
        
        % determine which cross best represents the transition based on type (i.e. high to low vs. low to high)
        slopeInt = 50; % interval for determining pos vs neg slope
        if isnan(leftCross) || leftCross <= slopeInt
            transCross = rightCross;
            leftEnd = leftCross;
            rightEnd = findCross(f1trace, f2trace, transCross, 'right');
        elseif isnan(rightCross) || rightCross >= length(f1trace) - slopeInt
            transCross = leftCross;
            rightEnd = rightCross;
            [leftEnd, ~] = findCross(f1trace, f2trace, transCross, 'left');
        else
            % check in higher frequency (generally has larger dynamic range)
            f2leftDiff = f2trace(leftCross - slopeInt) - f2trace(leftCross + slopeInt);
            f2rightDiff = f2trace(rightCross - slopeInt) - f2trace(rightCross + slopeInt);
            
            if tType == 5
                if f2leftDiff <= f2rightDiff
                    transCross = leftCross;
                    rightEnd = rightCross;
                    [leftEnd, ~] = findCross(f1trace, f2trace, transCross, 'left');
                else
                    transCross = rightCross;
                    leftEnd = leftCross;
                    rightEnd = findCross(f1trace, f2trace, transCross, 'right');
                end
            elseif tType == 6
                if f2leftDiff >= f2rightDiff
                    transCross = leftCross;
                    rightEnd = rightCross;
                    [leftEnd, ~] = findCross(f1trace, f2trace, transCross, 'left');
                else
                    transCross = rightCross;
                    leftEnd = leftCross;
                    rightEnd = findCross(f1trace, f2trace, transCross, 'right');
                end
            end
        end
        
        % no cross on left
        if isnan(leftEnd)
            leftEnd = 1;
        end
        
        % no cross on right
        if isnan(rightEnd)
            rightEnd = length(f1trace);
        end
        
        if tType == 5
            [~, f1PeakI] = nanmax(f1trace(leftEnd:transCross));
            f1PeakI = f1PeakI - 1 + leftEnd;
            [~, f1TroughI] = nanmin(f1trace(transCross:rightEnd));
            f1TroughI = f1TroughI - 1 + transCross;
            [~, f2PeakI] = nanmax(f2trace(transCross:rightEnd));
            f2PeakI = f2PeakI - 1 + transCross;
            [~, f2TroughI] = nanmin(f2trace(leftEnd:transCross));
            f2TroughI = f2TroughI - 1 + leftEnd;
            
        else
            [~, f1PeakI] = nanmax(f1trace(transCross:rightEnd));
            f1PeakI = f1PeakI - 1 + transCross;
            [~, f1TroughI] = nanmin(f1trace(leftEnd:transCross));
            f1TroughI = f1TroughI - 1 + leftEnd;
            [~, f2PeakI] = nanmax(f2trace(leftEnd:transCross));
            f2PeakI = f2PeakI - 1 + leftEnd;
            [~, f2TroughI] = nanmin(f2trace(transCross:rightEnd));
            f2TroughI = f2TroughI - 1 + transCross;
        end
        %                 if strcmp(parName, 'stratus132')
        %                                     figure;
        %                                     plot(f1trace)
        %                                     hold on
        %                                     plot(f2trace)
        %                                     plot(f1PeakI, f1trace(f1PeakI), 'ro');
        %                                     plot(f1TroughI, f1trace(f1TroughI), 'ro');
        %                                     plot(f2PeakI, f2trace(f2PeakI), 'bo');
        %                                     plot(f2TroughI, f2trace(f2TroughI), 'bo');
        %                                     hold off
        %                 end
        
        f1simAmps = [f1simAmps f1trace(f1PeakI) - f1trace(f1TroughI)];
        f2simAmps = [f2simAmps f2trace(f2PeakI) - f2trace(f2TroughI)];
        % error propagation
        f1PTerror = [f1PTerror sqrt((f1error(f1PeakI))^2 + (f1error(f1TroughI))^2)];
        f2PTerror = [f2PTerror sqrt((f2error(f2PeakI))^2 + (f2error(f2TroughI))^2)];
    end
end

% average peak to trough sim amplitude for 'high to low' and 'low
% to high' transitions
f1simAmp = nanmean(f1simAmps);
f2simAmp = nanmean(f2simAmps);
% error propagation
f1simAmp_error = sqrt((f1PTerror(1)/2)^2 + (f1PTerror(2)/2)^2);
f2simAmp_error = sqrt((f2PTerror(1)/2)^2 + (f2PTerror(2)/2)^2);

ntList = [];
ntLists = {ntList ntList};
normError = [];
normError = {normError normError}; % space to store standard error
s2dtrace = [];
d2strace = [];
if ~isnan(f1simAmp) % ensure that there are sim trials
    currTf1 = [];
    currTf2 = [];
    for tType = 1:2 % rivalry l2h, rivalry h2l
        if size(tLists{tType}.f1,1) > minInstances
            %for t = 1: size(tLists{tType}.f1,1)
            
            % individual transitions
            %                     currTf1 = tLists{tType}.f1(t,:);
            %                     currTf2 = tLists{tType}.f2(t,:);
            
            % averaged transition
            currTf1 = nanmean(tLists{tType}.f1,1);
            currTf2 = nanmean(tLists{tType}.f2,1);
            Tf1error = ste(tLists{tType}.f1);
            Tf2error = ste(tLists{tType}.f2);
            
            %normalize rivalry peaks using sim amplitudes
            
            % individual transitions
            %ntLists{tType}.f1(t,:) = (currTf1 - min(currTf1)) / f1simAmp;
            %ntLists{tType}.f2(t,:) = (currTf2 - min(currTf2)) / f2simAmp;
            
            % normalize the transition
            [minVal1, minI1] = min(currTf1);
            [minVal2, minI2] = min(currTf2);
            normedF1 = (currTf1 - minVal1) / f1simAmp;
            normedF2 = (currTf2 - minVal2) / f2simAmp;
            % error propagation
            normedF1_error = sqrt((1/f1simAmp)^2*((Tf1error).^2 + (Tf1error(minI1))^2) + ((currTf1 - minVal1)/f1simAmp^2*f1simAmp_error).^2);
            normedF2_error = sqrt((1/f2simAmp)^2*((Tf2error).^2 + (Tf2error(minI2))^2) + ((currTf2 - minVal2)/f2simAmp^2*f2simAmp_error).^2);
            
            %                     ntLists{tType}.f1 = (currTf1 - min(currTf1)) / f1simAmp;
            %                     ntLists{tType}.f2 = (currTf2 - min(currTf2)) / f2simAmp;
            
            % demean after normalizing
            ntLists{tType}.f1 = normedF1 - nanmean(normedF1);
            ntLists{tType}.f2 = normedF2 - nanmean(normedF2);
            % error propagation (final)
            normError{tType}.f1 = sqrt(normedF1_error.^2 + (1/length(normedF1))^2*sum(normedF1_error.^2));
            normError{tType}.f2 = sqrt(normedF2_error.^2 + (1/length(normedF2))^2*sum(normedF2_error.^2));
            
            
            %                     figure
            %                     plot(currTf1,'b--')
            %                     hold on
            %                     plot(ntLists{tType}.f1,'b')
            %                     plot(currTf2,'r--')
            %                     plot(ntLists{tType}.f2,'r')
            %                     title(parName)
            
            %end
            
            %             if ~isnan(ntLists{tType}.f1)
            %                 if tType == 1 % low to high
            %                     s2dtrace = [s2dtrace; ntLists{tType}.f2];
            %                     d2strace = [d2strace; ntLists{tType}.f1];
            %                 elseif tType ==2 % high to low
            %                     s2dtrace = [s2dtrace; ntLists{tType}.f1];
            %                     d2strace = [d2strace; ntLists{tType}.f2];
            %                 end
            %             end
            %         else
            %             if isempty(instExcluded) || ~strcmp(instExcluded(end), parName)
            %                 instExcluded = [instExcluded parName];
            %             end
        end
    end
end


%% Save collected transitions and peaks

% transitions
tLists = {l2hRival h2lRival l2lRival h2hRival l2hSim h2lSim l2hMarzRival h2lMarzRival l2lMarzRival h2hMarzRival};
tLabels = {'dartboard rivalry low to high' 'dartboard rivalry high to low' ...
    'dartboard rivalry low to low' 'dartboard rivalry high to high' ...
    'sim low to high' 'sim high to low' 'marz rivalry low to high' 'marz rivalry high to low' ...
    'marz rivalry low to low' 'marz rivalry high to high'};

save(['transitions/finalTransitions/' parName '_' normType '_' analFreqLabel '_' electrodeSet.name], 'tLists');

% normalized transitions
save(['transitions/finalTransitions_norm/' parName '_' normType '_' analFreqLabel '_' electrodeSet.name], 'ntLists');

% peaks
pLists = {slowPkRival fastPkRival slowPkSim fastPkSim slowPkMarzRival fastPkMarzRival};
pLabels = {'dartboard rivalry low dom' 'dartboard rivalry high dom' 'sim low dom' 'sim high dom' ...
    'marz rivalry low dom' 'marz rivalry high dom'};
pFileNames = {'dartRival_lowDom' 'dartRival_highDom' 'simRival_lowDom' 'simRival_highDom' ...
    'marzRival_lowDom' 'marzRival_highDom'};

save(['peaks/finalPeaks/' parName '_' normType], 'pLists');

%% Plot averaged transitions

if strcmp(parPlotOrNot, 'yes')
    %plot averaged transitions (low to high, high to low, low to low, high to high
    %for tType = 1:length(tLists)
    %     for tType = 1:2
    %         % must have a minimum number of transitions for each run type
    %         if size(tLists{tType}.f1, 1) < minInstances
    %             disp([parName 'does not have enough ' tLabels(tType) ' instances'])
    %             continue;
    %         end
    %
    %         meanF1trace = nanmean(tLists{tType}.f1,1);
    %         errorF1trace = ste(tLists{tType}.f1);
    %
    %         meanF2trace = nanmean(tLists{tType}.f2,1);
    %         errorF2trace = ste(tLists{tType}.f2);
    %
    %         figure
    %         title([parName ' ' tLabels{tType} ': ' num2str(size(tLists{tType}.f1, 1)) ' transitions averaged'])
    %         hold on
    %         mseb([(-transHalf):(1/512):transHalf],[meanF1trace; meanF2trace],[errorF1trace; errorF2trace],[],1);
    %         ylim([-1.5 1.5])
    %         vline(0,'k','transition')
    %         xlabel('Time from button press (s)');
    %         ylabel('Amplitude');
    %         legend('Low Frequency', 'High Frequency');
    %     end
    
    % dartRival
    % must have a minimum number of transitions
    if size([tLists{1}.f2;tLists{2}.f1]) >= minInstances
        
        % plot dartRival transitions: l2h and h2l averaged
        s2dTrace = nanmean([tLists{1}.f2; tLists{2}.f1],1);
        s2dError = ste([tLists{1}.f2; tLists{2}.f1]);
        d2sTrace = nanmean([tLists{1}.f1; tLists{2}.f2],1);
        d2sError = ste([tLists{1}.f1; tLists{2}.f2]);
        
        figure
        title(['       ' parName ' dartRival: ' num2str(size([tLists{1}.f1; tLists{2}.f1], 1)) ' transitions     AB'], 'FontSize', 16)
        hold on
        mseb([(-transHalf):(1/sampRate):transHalf],[s2dTrace; d2sTrace],[s2dError; d2sError],[],1);
        ylim([-0.8 0.8])
        vline(0,'k','transition')
        xlabel('Time from button press (s)');
        ylabel('Amplitude');
        legend('Frequency 1', 'Frequency 2');
        
        % save dartRival transitions figure
        figName = [transPlotDir 'Individual\' parName '_dartRival'];
        saveas(gcf,figName,'jpg')
    else
        disp([parName ' does not have enough transitions'])
    end
    
    % dartSim
    % plot dartSim transitions: l2h and h2l averaged
    s2dTrace = nanmean([tLists{5}.f2; tLists{6}.f1],1);
    s2dError = ste([tLists{5}.f2; tLists{6}.f1]);
    d2sTrace = nanmean([tLists{5}.f1; tLists{6}.f2],1);
    d2sError = ste([tLists{5}.f1; tLists{6}.f2]);
    
    figure
    title(['            ' parName ' dartSim: ' num2str(size([tLists{5}.f1; tLists{6}.f1], 1)) ' transitions          AB'], 'FontSize', 16)
    hold on
    mseb([(-transHalf):(1/sampRate):transHalf],[s2dTrace; d2sTrace],[s2dError; d2sError],[],1);
    ylim([-0.8 0.8])
    vline(0,'k','transition')
    xlabel('Time from button press (s)');
    ylabel('Amplitude');
    legend('Frequency 1', 'Frequency 2');
    
    % save dartSim transitions figure
    figName = [transPlotDir 'Individual\' parName '_dartSim'];
    saveas(gcf,figName,'jpg')
    
    % dartRival normalized to dartSim
    if ~isempty(ntLists{1}) || ~isempty(ntLists{2})
        
        % plot dartRival transitions: l2h and h2l averaged
        s2dTrace = nanmean([ntLists{1}.f2; ntLists{2}.f1],1);
        s2dError = sqrt((normError{1}.f2/2).^2 + (normError{2}.f1/2).^2);
        %s2dError = ste([ntLists{1}.f2; ntLists{2}.f1]);
        d2sTrace = nanmean([ntLists{1}.f1; ntLists{2}.f2],1);
        d2sError = sqrt((normError{1}.f1/2).^2 + (normError{2}.f2/2).^2);
        %d2sError = ste([ntLists{1}.f1; ntLists{2}.f2]);
        
        figure
        title([' Normalized Average Transition: ' parName '     .'], 'FontSize', 16)
        hold on
        mseb([(-transHalf):(1/sampRate):transHalf],[s2dTrace; d2sTrace],[s2dError; d2sError],[],1);
        %mseb([(-transHalf):(1/sampRate):transHalf],[s2dTrace; d2sTrace],[normError1; normError2],[],1);
        ylim([-1.9 1.9])
        vline(0,'k','transition')
        xlabel('Time from button press (s)');
        ylabel('Amplitude');
        legend('Frequency 1', 'Frequency 2');
        legend boxoff
        
        % save dartRival transitions figure
        figName = [transPlotDir 'Individual\' parName '_dartRival_norm'];
        saveas(gcf,figName,'jpg')
    else
        disp([parName 'does not have enough transitions (for normalized transitions)'])
    end
    
    %% plot averaged peaks
    
    %     for pType = 1:length(pLists) % 1:2
    %
    %         % must have a minimum number of peaks for each run type
    %         if size(pLists{pType}.f1, 1) < minInstances
    %             disp([parName 'does not have enough ' pLabels(pType) ' peaks'])
    %             continue;
    %         end
    %
    %         meanF1trace = nanmean(pLists{pType}.f1,1);
    %         errorF1trace = ste(pLists{pType}.f1);
    %
    %         meanF2trace = nanmean(pLists{pType}.f2,1);
    %         errorF2trace = ste(pLists{pType}.f2);
    %
    %         figure
    %         title([parName ' ' pLabels{pType} ': ' num2str(size(pLists{pType}.f1, 1)) ' peaks'], 'FontSize', 14)
    %         hold on
    %         mseb([(-peakHalf):(1/512):peakHalf],[meanF1trace; meanF2trace],[errorF1trace; errorF2trace],[],1);
    %         vline(0,'k','peak')
    %         xlabel('Time from peak (s)');
    %         ylabel('Amplitude');
    %         legend('Low Frequency', 'High Frequency');
    %     end
end

end


