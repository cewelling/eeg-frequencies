function makeGroupPlot()
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

clearvars

% load parameters
groupAnalysisParams

numFormat = '%02d'; % for participant number

tLabels = {'dartboard rivalry low to high' 'dartboard rivalry high to low' 'dartboard rivalry low to low' 'dartboard rivalry high to high' 'sim low to high' 'sim high to low' ...
    'marz rivalry low to high' 'marz rivalry high to low' 'marz rivalry low to low' 'marz rivalry high to high' 'norm rivalry low to high' 'norm rivalry high to low' 'dartboard rivalry' 'dartboard sim' 'normalized transition'};
tFileNames = {'dartRival_l2h' 'dartRival_h2l' 'dartRival_l2l' 'dartRival_h2h' 'dartSim_l2h' 'dartSim_h2l' 'marzRival_l2h' 'marzRival_h2l' 'marzRival_l2l' 'marzRival_h2h' 'normRival_l2h' 'normRival_h2l' 'dartRival' 'dartSim' 'dartNorm'};

pLabels = {'dartboard rivalry low dom' 'dartboard rivalry high dom' 'sim low dom' 'sim high dom' ...
    'marz rivalry low dom' 'marz rivalry high dom' 'norm rivalry low dom' 'norm rivalry high dom' 'dartboard rivalry' 'dartboard sim' 'normalized peak'};
pFileNames = {'dartRival_lowDom' 'dartRival_highDom' 'dartSim_lowDom' 'dartSim_highDom' ...
    'marzRival_lowDom' 'marzRival_highDom' 'normRival_lowDom' 'normRival_highDom' 'dartRival' 'dartSim' 'dartNorm'};

yesOrNo = input('Use last peaks and transitions stored?', 's');
if strcmp(yesOrNo, 'yes')
    loaded = 1;
else
    loaded = 0;
end

for iGroup = 1:length(groupParNums)
    
    % allocate space to store transitions
    l2hRival.f1 = []; l2hRival.f2 = [];
    h2lRival.f1 = []; h2lRival.f2 = [];
    l2lRival.f1 = []; l2lRival.f2 = [];
    h2hRival.f1 = []; h2hRival.f2 = [];
    l2hSim.f1 = []; l2hSim.f2 = [];
    h2lSim.f1 = []; h2lSim.f2 = [];
    l2hMarzRival.f1 = []; l2hMarzRival.f2 = [];
    h2lMarzRival.f1 = []; h2lMarzRival.f2 = [];
    l2lMarzRival.f1 = []; l2lMarzRival.f2 = [];
    h2hMarzRival.f1 = []; h2hMarzRival.f2 = [];
    l2hNorm.f1 = []; l2hNorm.f2 = [];
    h2lNorm.f1 = []; h2lNorm.f2 = [];
    dartRival.f1 = []; dartRival.f2 = [];
    dartSim.f1 = []; dartSim.f2 = [];
    dartNorm.f1 = []; dartNorm.f2 = [];
    
    groupTLists = {l2hRival h2lRival l2lRival h2hRival l2hSim h2lSim l2hMarzRival h2lMarzRival l2lMarzRival h2hMarzRival l2hNorm h2lNorm dartRival dartSim dartNorm};
    
    % allocate space to store peaks
    slowPkRival.f1 = []; slowPkRival.f2 = [];
    fastPkRival.f1 = []; fastPkRival.f2 = [];
    slowPkSim.f1 = []; slowPkSim.f2 = [];
    fastPkSim.f1 = []; fastPkSim.f2 = [];
    slowPkMarzRival.f1 = []; slowPkMarzRival.f2 = [];
    fastPkMarzRival.f1 = []; fastPkMarzRival.f2 = [];
    slowPkNorm.f1 = []; slowPkNorm.f2 = [];
    fastPkNorm.f1 = []; fastPkNorm.f2 = [];
    dartRival.f1 = []; dartRival.f2 = [];
    dartSim.f1 = []; dartSim.f2 = [];
    dartNorm.f1 = []; dartNorm.f2 = [];
    
    groupPLists = {slowPkRival fastPkRival slowPkSim fastPkSim slowPkMarzRival fastPkMarzRival slowPkNorm fastPkNorm dartRival dartSim dartNorm};
    
    % Store sim amps and rivalry maxima to compare
    f1simAmpM = [];
    f2simAmpM = [];
    f1rivMaxes = [];
    f2rivMaxes = [];
    
    % Store participants excluded due to low number of transitions
    instExcluded = {};
    
    for parNum = groupParNums{iGroup};
        parName = [groupCodes{iGroup} num2str(parNum, numFormat)];
        
        % Get EEG file info
        EEGfiles = dir([eegDir parName '*']);
        % ...if the participant exists...
        if isempty(EEGfiles)
            continue;
        end
        date = strtok(EEGfiles(1).date);
        date = datestr(date, dateFormat);
        
        % get participant's transitions and peaks
        if loaded
            try
                if altPeak
                    load(['peaks/peaks2/' parName '_' normType '.mat'])
                else
                    load(['peaks/finalPeaks/' parName '_' normType '.mat'])
                end
                load(['transitions/finalTransitions/' parName '_' normType '.mat'])
            catch
                [ tLists, ~, pLists, ~] = makeParPlot(parName, date);
            end
        else
            [ tLists, ~, pLists, ~] = makeParPlot(parName, date);
        end
        
        %% Collect dartRival transitions normalized to dartSim amplitudes
        
        % Find average sim amplitude
        
        %collect sim maxes and mins
        %         f1maxes = [];
        %         f1mins = [];
        %         f2maxes = [];
        %         f2mins = [];
        
        f1simAmps = [];
        f2simAmps = [];
        for tType = 5:6 % 5: sim low to high; 6: sim high to low
            
            % average transitions
            f1trace = nanmean(tLists{tType}.f1,1);
            f2trace = nanmean(tLists{tType}.f2,1);
            
            if ~isnan(f1trace) % no f1 transitions means no f2 transitions either
                
                % find peak and trough of averaged transition
                
                % demean time courses
                f1trace = f1trace - nanmean(f1trace);
                f2trace = f2trace - nanmean(f2trace);
                
                % identify where they cross on both sides of button press
                [leftCross, leftDist] = findCross(f1trace, f2trace, buttonPress, 'left');
                [rightCross, rightDist] = findCross(f1trace, f2trace, buttonPress, 'right');
                
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
                    try
                    [~, f1PeakI] = nanmax(f1trace(leftEnd:transCross));
                    catch
                        disp('hi');
                    end
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
                
                %                     figure;
                %                     plot(f1trace)
                %                     hold on
                %                     plot(f2trace)
                %                     plot(f1PeakI, f1trace(f1PeakI), 'ro');
                %                     plot(f1TroughI, f1trace(f1TroughI), 'ro');
                %                     plot(f2PeakI, f2trace(f2PeakI), 'bo');
                %                     plot(f2TroughI, f2trace(f2TroughI), 'bo');
                %                     hold off
                
                f1simAmps = [f1simAmps f1trace(f1PeakI) - f1trace(f1TroughI)];
                f2simAmps = [f2simAmps f2trace(f2PeakI) - f2trace(f2TroughI)];
                %             f1maxes = [f1maxes; max(tLists{tType}.f1, [],2)];
                %             f1mins = [f1mins; min(tLists{tType}.f1, [], 2)];
                %             f2maxes = [f2maxes; max(tLists{tType}.f2, [],2)];
                %             f2mins = [f2mins; min(tLists{tType}.f2, [], 2)];
            end
        end
        
        % calculate average peak to trough sim amplitude
        f1simAmp = nanmean(f1simAmps);
        f2simAmp = nanmean(f2simAmps);
        %         f1simAmp = nanmean(f1maxes) - nanmean(f1mins);
        %         f2simAmp = nanmean(f2maxes) - nanmean(f2mins);
        
        ntList = [];
        ntLists = {ntList ntList};
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
                    
                    %normalize rivalry peaks using sim amplitudes
                    
                    % individual transitions
                    %ntLists{tType}.f1(t,:) = (currTf1 - min(currTf1)) / f1simAmp;
                    %ntLists{tType}.f2(t,:) = (currTf2 - min(currTf2)) / f2simAmp;
                    
                    % averaged transition
                    ntLists{tType}.f1 = (currTf1 - min(currTf1)) / f1simAmp;
                    ntLists{tType}.f2 = (currTf2 - min(currTf2)) / f2simAmp;
                    
                    f1simAmpM(parNum) = f1simAmp;
                    f2simAmpM(parNum) = f2simAmp;
                    f1rivMaxes(parNum) = max(currTf1) - min(currTf1);
                    f2rivMaxes(parNum) = max(currTf2) - min(currTf2);
                    %end
                    
                    if ~isnan(ntLists{tType}.f1)
                        if tType == 1 % low to high
                            s2dtrace = [s2dtrace; ntLists{tType}.f2];
                            d2strace = [d2strace; ntLists{tType}.f1];
                        elseif tType ==2 % high to low
                            try
                                s2dtrace = [s2dtrace; ntLists{tType}.f1];
                            catch
                                disp('hi')
                            end
                            d2strace = [d2strace; ntLists{tType}.f2];
                        end
                    end
                else
                    if isempty(instExcluded) || ~strcmp(instExcluded(end), parName)
                        instExcluded = [instExcluded parName];
%                     elseif ~strcmp(instExcluded(end), parName)
%                         instExcluded = [instExcluded parName];
                    end
                end
            end
        end
        
        % align averaged normalized transitions to transition CROSS
        meanS2D = nanmean(s2dtrace,1);
        meanD2S = nanmean(d2strace,1);
        
        % identify where traces cross on both sides of button press
        [leftCross, leftDist] = findCross(meanS2D, meanD2S, buttonPress, 'left');
        [rightCross, rightDist] = findCross(meanS2D, meanD2S, buttonPress, 'right');
        
        % determine which cross best represents the transition based on type (i.e. high to low vs. low to high)
        slopeInt = 50; % interval for determining pos vs neg slope
        if isnan(leftCross) || leftCross <= slopeInt
            transCross = rightCross;
            dist = rightDist;
        elseif isnan(rightCross) || rightCross >= length(f1trace) - slopeInt
            transCross = leftCross;
            dist = -leftDist;
        else
            % check in just one frequency
            leftDiff = meanS2D(leftCross - slopeInt) - meanS2D(leftCross + slopeInt);
            rightDiff = meanS2D(rightCross - slopeInt) - meanS2D(rightCross + slopeInt);
            
            if leftDiff <= rightDiff
                transCross = leftCross;
                dist = -leftDist;
            else
                transCross = rightCross;
                dist = rightDist;
            end
        end
        
        % align to transition cross, pad with NaNs, store
        if ~isnan(meanS2D)
            alignedS2D = [nan(1,sampRate * 5 - dist), meanS2D, nan(1,sampRate * 5 + dist)];
            groupTLists{15}.f1 = [groupTLists{15}.f1; alignedS2D];
        end
        if ~isnan(meanD2S)
            alignedD2S = [nan(1,sampRate * 5 - dist), meanD2S, nan(1,sampRate * 5 + dist)];
            groupTLists{15}.f2 = [groupTLists{15}.f2; alignedD2S];
        end
        
        %         if ~isnan(s2dtrace)
        %             groupTLists{15}.f1 = [groupTLists{15}.f1; nanmean(s2dtrace,1)];
        %         end
        %         if ~isnan(d2strace)
        %             groupTLists{15}.f2 = [groupTLists{15}.f2; nanmean(d2strace,1)];
%         end
        
        save(['transitions/finalTransitions_norm/' parName '_' normType], 'ntLists');
        
        if ~isempty(ntLists)
            % mean of participant's normalized transitions
            if ~isempty(ntLists{1})
                l2hRivf1 = nanmean(ntLists{1}.f1,1);
                l2hRivf2 = nanmean(ntLists{1}.f2,1);
            else
                l2hRivf1 = [];
                l2hRivf2 = [];
            end
            if ~isempty(ntLists{2})
                h2lRivf1 = nanmean(ntLists{2}.f1,1);
                h2lRivf2 = nanmean(ntLists{2}.f2,1);
            else
                h2lRivf1 = [];
                h2lRivf2 = [];
            end
            
            % collect normalized transitions
            if size(tLists{1}.f1,1) > minInstances
                groupTLists{11}.f1 = [groupTLists{11}.f1; l2hRivf1];
                groupTLists{11}.f2 = [groupTLists{11}.f2; l2hRivf2];
            end
            if size(tLists{2}.f1,1) > minInstances
                groupTLists{12}.f1 = [groupTLists{12}.f1; h2lRivf1];
                groupTLists{12}.f2 = [groupTLists{12}.f2; h2lRivf2];
            end
        end
        
        
        %% Collect transitions
        s2dtraceR = [];
        d2straceR = [];
        s2dtraceS = [];
        d2straceS = [];
        for tType = 1:length(tLists) - 2 % last ones are for h2l and l2h averaged
            if size(tLists{tType}.f1,1) > minInstances
                
                % mean of participant's transitions
                f1mean = nanmean(tLists{tType}.f1,1);
                f2mean = nanmean(tLists{tType}.f2,1);
                if isnan(f1mean)
                    f1mean = [];
                end
                if isnan(f2mean)
                    f2mean = [];
                end
                
                % store transitions
                groupTLists{tType}.f1 = [groupTLists{tType}.f1; f1mean];
                groupTLists{tType}.f2 = [groupTLists{tType}.f2; f2mean];
                
                if tType == 1  % rivalry 
                    % supp to dom trace
                    s2dtraceR = [s2dtraceR; f2mean];
                    % dom to supp trace
                    d2straceR = [d2straceR; f1mean];
                elseif tType == 2 % rivalry
                    % supp to dom trace
                    s2dtraceR = [s2dtraceR; f1mean];
                    % dom to supp trace
                    d2straceR = [d2straceR; f2mean];
                elseif tType == 5 % sim
                    % supp to dom trace
                    s2dtraceS = [s2dtraceS; f2mean];
                    % dom to supp trace
                    d2straceS = [d2straceS; f1mean];
                elseif tType == 6 % sim
                    % supp to dom trace
                    s2dtraceS = [s2dtraceS; f1mean];
                    % dom to supp trace
                    d2straceS = [d2straceS; f2mean];
                end
            end
        end
        
        % average low to high and high to low

        if ~isnan(s2dtraceR)
            groupTLists{13}.f1 = [groupTLists{13}.f1; nanmean(s2dtraceR,1)];
        end
        if ~isnan(d2straceR)
            groupTLists{13}.f2 = [groupTLists{13}.f2; nanmean(d2straceR,1)];
        end
        if ~isnan(s2dtraceS)
            groupTLists{14}.f1 = [groupTLists{14}.f1; nanmean(s2dtraceS,1)];
        end
        if ~isnan(d2straceS)
            groupTLists{14}.f2 = [groupTLists{14}.f2; nanmean(d2straceS,1)];
        end
        
        %% Collect dartRival peaks normalized to dartSim amplitudes
        
        % Find average sim amplitude
        
        f1simAmps = [];
        f2simAmps = [];
        
        % collect sim maxes and minx
        for pType = 3:4 % 3: sim low to high; 4: sim high to low
            
            % average transitions
            f1trace = nanmean(pLists{pType}.f1);
            f2trace = nanmean(pLists{pType}.f2);
            
            if ~isnan(f1trace) % no f1 transitions means no f2 transitions either
                
                % find peak and trough of averaged transition
                
                if pType == 3 % 3: sim low dom
                    f1lowMax = max(f1trace);
                    f1lowMin = min(f1trace);
                    f2lowMax = max(f2trace);
                    f2lowMin = min(f2trace);
                elseif pType == 4 % 4: sim high dom
                    f1highMax = max(f1trace);
                    f1highMin = min(f1trace);
                    f2highMax = max(f2trace);
                    f2highMin = min(f2trace);
                end
                
                
                %             if pType == 3 % 3: sim low dom
                %                 f1lowMaxes = max(pLists{pType}.f1, [],2);
                %                 f1lowMins = min(pLists{pType}.f1, [], 2);
                %                 f2lowMaxes = max(pLists{pType}.f2, [],2);
                %                 f2lowMins = min(pLists{pType}.f2, [], 2);
                %             elseif pType == 4 % 4: sim high dom
                %                 f1highMaxes = max(pLists{pType}.f1, [],2);
                %                 f1highMins = min(pLists{pType}.f1, [], 2);
                %                 f2highMaxes = max(pLists{pType}.f2, [],2);
                %                 f2highMins = min(pLists{pType}.f2, [], 2);
            end
        end
        
        % calculate average peak to trough sim amplitude
        f1lsimAmp = f1lowMax - f1lowMin;
        f1hsimAmp = f1highMax - f1highMin;
        f2lsimAmp = f2lowMax - f2lowMin;
        f2hsimAmp = f2highMax - f2highMin;
        
%         % calculate average peak to trough sim amplitude
%         f1lsimAmp = nanmean(f1lowMaxes) - nanmean(f1lowMins);
%         f1hsimAmp = nanmean(f1highMaxes) - nanmean(f1highMins);
%         f2lsimAmp = nanmean(f2lowMaxes) - nanmean(f2lowMins);
%         f2hsimAmp = nanmean(f2highMaxes) - nanmean(f2highMins);
        
        npList = [];
        npLists = {npList npList};
        domtrace = [];
        supptrace = [];
        if ~isnan(f1lsimAmp) % ensure that there are sim trials
            currPeakf1 = [];
            currPeakf2 = [];
            for pType = 1:2 % rivalry low dom, rivalry high dom
                %                 for p = 1: size(pLists{pType}.f1,1)
                
                % first average participant's peaks
                currPeakf1 = nanmean(pLists{pType}.f1);
                currPeakf2 = nanmean(pLists{pType}.f2);
                
                %                     currPeakf1 = pLists{pType}.f1(p,:);
                %                     currPeakf2 = pLists{pType}.f2(p,:);
                
                %normalize rivalry peaks using sim amplitudes
                    if pType == 1
                        npLists{pType}.f1 = (currPeakf1 - min(currPeakf1)) / f1lsimAmp;
                        npLists{pType}.f2 = (currPeakf2 - min(currPeakf2)) / f2lsimAmp;
                    else
                        npLists{pType}.f1 = (currPeakf1 - min(currPeakf1)) / f1hsimAmp;
                        npLists{pType}.f2 = (currPeakf2 - min(currPeakf2)) / f2hsimAmp;
                    end
                    %                     if pType == 1
                    %                         npLists{pType}.f1(p,:) = (currPeakf1 - min(currPeakf1)) / f1lsimAmp;
                    %                         npLists{pType}.f2(p,:) = (currPeakf2 - min(currPeakf2)) / f2lsimAmp;
                    %                     else
                    %                         npLists{pType}.f1(p,:) = (currPeakf1 - min(currPeakf1)) / f1hsimAmp;
                    %                         npLists{pType}.f2(p,:) = (currPeakf2 - min(currPeakf2)) / f2hsimAmp;
                    %                     end
                    %                 end
                    
                    if pType == 1 % low to high
                        domtrace = [domtrace; npLists{pType}.f1];
                        supptrace = [supptrace; npLists{pType}.f2];
                    elseif tType ==2 % high to low
                        domtrace = [domtrace; npLists{pType}.f2];
                        supptrace = [supptrace; npLists{pType}.f1];
                    end                    
            end
        end
        
        % average high dom and low dom peaks
        if ~isnan(domtrace)
            groupPLists{11}.f1 = [groupPLists{11}.f1; nanmean(domtrace,1)];
        end
        if ~isnan(supptrace)
            groupPLists{11}.f2 = [groupPLists{11}.f2; nanmean(supptrace,1)];
        end
        
        save(['peaks/finalPeaks_norm/' parName '_' normType], 'npLists');
        
        if ~isempty(npLists) % snr too low
            
            % mean of participant's normalized peaks
            if ~isempty(npLists{1})
                ldomRivf1 = nanmean(npLists{1}.f1,1);
                ldomRivf2 = nanmean(npLists{1}.f2,1);
            else
                ldomRivf1 = [];
                ldomRivf2 = [];
            end
            if ~isempty(npLists{2})
                hdomRivf1 = nanmean(npLists{2}.f1,1);
                hdomRivf2 = nanmean(npLists{2}.f2,1);
            else
                hdomRivf1 = [];
                hdomRivf2 = [];
            end
            
            % collect peaks
            if size(pLists{1}.f1,1) > minInstances
                groupPLists{7}.f1 = [groupPLists{7}.f1; ldomRivf1];
                groupPLists{7}.f2 = [groupPLists{7}.f2; ldomRivf2];
            end
            if size(pLists{2}.f1,1) > minInstances
                groupPLists{8}.f1 = [groupPLists{8}.f1; hdomRivf1];
                groupPLists{8}.f2 = [groupPLists{8}.f2; hdomRivf2];
            end
        end
        %normAmp = mean([max(l2hRivf1) max(h2lRivf1) max(l2hRivf2) max(h2lRivf2)]);
        %save([indicesDir 'normAmps/' parName], 'normAmp')
        
        %% Collect peaks
        
        domtraceR = [];
        domtraceS = [];
        supptraceR = [];
        supptraceS = [];
        for pType = 1:length(pLists)
            
            if size(pLists{pType}.f1,1) > minInstances
                % mean of participant's transitions
                f1mean = nanmean(pLists{pType}.f1,1);
                f2mean = nanmean(pLists{pType}.f2,1);
                %                 if isnan(f1mean)
                %                     f1mean = [];
                %                 end
                %                 if isnan(f2mean)
                %                     f2mean = [];
                %                 end
                
                % store peaks
                groupPLists{pType}.f1 = [groupPLists{pType}.f1; f1mean];
                groupPLists{pType}.f2 = [groupPLists{pType}.f2; f2mean];
                
                 if pType == 1  % rivalry 
                    domtraceR = [domtraceR; f1mean];
                    supptraceR = [supptraceR; f2mean];
                elseif pType == 2 % rivalry
                    domtraceR = [domtraceR; f2mean];
                    supptraceR = [supptraceR; f1mean];
                elseif pType == 3 % sim
                    domtraceS = [domtraceS; f1mean];
                    supptraceS = [supptraceS; f2mean];
                elseif pType == 4 % sim
                    domtraceS = [domtraceS; f2mean];
                    supptraceS = [supptraceS; f1mean];
                end
            end
        end
        
                % average high and low dom peaks
        if ~isnan(domtraceR)
            try
            groupPLists{9}.f1 = [groupPLists{9}.f1; nanmean(domtraceR,1)];
            catch
                disp('hi')
            end
        end
        if ~isnan(supptraceR)
            groupPLists{9}.f2 = [groupPLists{9}.f2; nanmean(supptraceR,1)];
        end
        if ~isnan(domtraceS)
            groupPLists{10}.f1 = [groupPLists{10}.f1; nanmean(domtraceS,1)];
        end
        if ~isnan(supptraceS)
            groupPLists{10}.f2 = [groupPLists{10}.f2; nanmean(supptraceS,1)];
        end
        
    end
    
    %% Plotting
    
    % transitions
    %for tType = 1:length(groupTLists)
    for tType = [13 15]; %[11:12 15] %11:12
        
        %         if tType == 13
        %             meanF1trace = nanmean([groupTLists{1}.f1; groupTLists{2}.f2],1);
        %             errorF1trace = ste([groupTLists{1}.f1; groupTLists{2}.f2]);
        %
        %             meanF2trace = nanmean([groupTLists{1}.f2; groupTLists{2}.f1],1);
        %             errorF2trace = ste([groupTLists{1}.f2; groupTLists{2}.f1]);
        %         else
        % must have a minimum number of transitions for each run type
        if size(groupTLists{tType}.f1, 1) < minPars
            disp([analGroupIDs{iGroup} 'does not have enough ' tLabels(tType) ' participants'])
            continue;
        end
        
        meanF1trace = nanmean(groupTLists{tType}.f1,1);
        errorF1trace = ste(groupTLists{tType}.f1);
        
        meanF2trace = nanmean(groupTLists{tType}.f2,1);
        errorF2trace = ste(groupTLists{tType}.f2);
        %         end
        
        figure
       %for i = 1:size(groupTLists{tType}.f1,1)
        plot(groupTLists{tType}.f1', 'b')
        %plot(groupTLists{tType}.f1(i,:), 'b')
        hold on
        plot(groupTLists{tType}.f2', 'r')
        %plot(groupTLists{tType}.f2(i,:), 'r')
        %end
        
        figure
        title([analGroupIDs{iGroup} ' ' tLabels{tType} ': ' num2str(size(groupTLists{tType}.f1, 1)) ' participants averaged'])
        hold on
        if tType == 15
            mseb([(-2*transHalf):(1/512):2*transHalf],[meanF1trace; meanF2trace],[errorF1trace; errorF2trace],[],1);
        else
            mseb([(-transHalf):(1/512):transHalf],[meanF1trace; meanF2trace],[errorF1trace; errorF2trace],[],1);
        end
        
        if tType == 11 || tType == 12 || tType == 15
            ylim([0 1.5])
        else
            ylim([-.4 .4]) %[-1 2])
        end
        vline(0,'k','transition')
        xlabel('Time from button press (s)');
        ylabel('Amplitude');
        
        if tType >= 13
            legend('Frequency 1', 'Frequency 2')
        else
            legend('Low Frequency', 'High Frequency');
        end
        
        if allTrans
            figName = [groupPlotDir analGroupIDs{iGroup} '_' tFileNames{tType} '_' num2str(size(groupTLists{tType}.f1,1)) 'pars_ALL.jpg'];
        else
            figName = [groupPlotDir analGroupIDs{iGroup} '_' tFileNames{tType} '_' num2str(size(groupTLists{tType}.f1,1)) 'pars.jpg'];
        end
        saveas(gcf,figName,'jpg')
    end
    
    % peaks
%         for pType = 1:length(groupPLists)
%             %for pType = 7:8
%     
%             % must have a minimum number of transitions for each run type
%             if size(groupPLists{pType}.f1, 1) < minInstances
%                 disp([analGroupIDs{iGroup} 'does not have enough ' pLabels(pType) ' instances'])
%                 continue;
%             end
%     
%             meanF1trace = nanmean(groupPLists{pType}.f1,1);
%             errorF1trace = ste(groupPLists{pType}.f1);
%     
%             meanF2trace = nanmean(groupPLists{pType}.f2,1);
%             errorF2trace = ste(groupPLists{pType}.f2);
%     
%             figure
%             title([analGroupIDs{iGroup} ' ' pLabels{pType} ': ' num2str(size(groupPLists{pType}.f1, 1)) ' participants averaged'])
%             hold on
%             mseb([(-peakHalf):(1/512):peakHalf],[meanF1trace; meanF2trace],[errorF1trace; errorF2trace],[],1);
%             if pType == 7 || pType == 8 || pType == 11
%                 ylim([0 1.5])
%             else
%                 ylim([-.4 .6]); %[-1 1.9])
%             end
%             vline(0,'k','peak')
%             xlabel('Time from peak (s)');
%             ylabel('Amplitude');
%             if pType >= 9
%                 legend('Frequency 1', 'Frequency 2');
%             else
%                 legend('Low Frequency', 'High Frequency');
%             end
%     
%             if altPeak
%                 figName = [groupPlotDir analGroupIDs{iGroup} '_' pFileNames{pType} 'Peaks_' num2str(size(groupPLists{pType}.f1,1)) 'pars_ALT.jpg'];
%             else
%                 figName = [groupPlotDir analGroupIDs{iGroup} '_' pFileNames{pType} 'Peaks_' num2str(size(groupPLists{pType}.f1,1)) 'pars.jpg'];
%             end
%             saveas(gcf,figName,'jpg')
%         end
end
end


