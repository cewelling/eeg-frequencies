 function makeGroupPlot()
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% load parameters
groupAnalysisParams

numFormat = '%02d'; % for participant number

tLabels = {'dartboard rivalry low to high' 'dartboard rivalry high to low' 'sim low to high' 'sim high to low' ...
    'marz rivalry low to high' 'marz rivalry high to low' 'norm rivalry low to high' 'norm rivalry high to low'};
tFileNames = {'dartRival_l2h' 'dartRival_h2l' 'dartSim_l2h' 'dartSim_h2l' 'marzRival_l2h' 'marzRival_h2l' 'normRival_l2h' 'normRival_h2l'};

pLabels = {'dartboard rivalry low dom' 'dartboard rivalry high dom' 'sim low dom' 'sim high dom' ...
    'marz rivalry low dom' 'marz rivalry high dom' 'norm rivalry low dom' 'norm rivalry high dom'};
pFileNames = {'dartRival_lowDom' 'dartRival_highDom' 'dartSim_lowDom' 'dartSim_highDom' ...
    'marzRival_lowDom' 'marzRival_highDom' 'normRival_lowDom' 'normRival_highDom'};

for iGroup = 1:length(groupParNums)
    
    % allocate space to store transitions
    l2hRival.f1 = []; l2hRival.f2 = [];
    h2lRival.f1 = []; h2lRival.f2 = [];
    l2hSim.f1 = []; l2hSim.f2 = [];
    h2lSim.f1 = []; h2lSim.f2 = [];
    l2hMarzRival.f1 = []; l2hMarzRival.f2 = [];
    h2lMarzRival.f1 = []; h2lMarzRival.f2 = [];
    l2hNorm.f1 = []; l2hNorm.f2 = [];
    h2lNorm.f1 = []; h2lNorm.f2 = [];
    groupTLists = {l2hRival h2lRival l2hSim h2lSim l2hMarzRival h2lMarzRival l2hNorm h2lNorm};
    
    % allocate space to store peaks
    slowPkRival.f1 = []; slowPkRival.f2 = [];
    fastPkRival.f1 = []; fastPkRival.f2 = [];
    slowPkSim.f1 = []; slowPkSim.f2 = [];
    fastPkSim.f1 = []; fastPkSim.f2 = [];
    slowPkMarzRival.f1 = []; slowPkMarzRival.f2 = [];
    fastPkMarzRival.f1 = []; fastPkMarzRival.f2 = [];
    slowPkNorm.f1 = []; slowPkNorm.f2 = [];
    fastPkNorm.f1 = []; fastPkNorm.f2 = [];
    groupPLists = {slowPkRival fastPkRival slowPkSim fastPkSim slowPkMarzRival fastPkMarzRival slowPkNorm fastPkNorm};
    
    for parNum = groupParNums{iGroup};
        parName = [groupCodes{iGroup} num2str(parNum, numFormat)];
        
        % load participant's transitions (if participant exists)
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
        
        % load participant's peaks (if participant exists)
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
        
        %% Collect dartRival transitions normalized to dartSim amplitudes
        
        % Find average sim amplitude
        
        %collect sim maxes and mins
        f1maxes = [];
        f1mins = [];
        f2maxes = [];
        f2mins = [];
        for tType = 3:4 % 3: sim low to high; 4: sim high to low
            f1maxes = [f1maxes; max(tLists{tType}.f1, [],2)];
            f1mins = [f1mins; min(tLists{tType}.f1, [], 2)];
            f2maxes = [f2maxes; max(tLists{tType}.f2, [],2)];
            f2mins = [f2mins; min(tLists{tType}.f2, [], 2)];
        end
        
        % calculate average peak to trough sim amplitude
        f1simAmp = nanmean(f1maxes) - nanmean(f1mins);
        f2simAmp = nanmean(f2maxes) - nanmean(f2mins);
        
        % mean of participant's rivalry transitions
        l2hRivf1 = nanmean(tLists{1}.f1,1);
        h2lRivf1 = nanmean(tLists{2}.f1,1);
        l2hRivf2 = nanmean(tLists{1}.f2,1);
        h2lRivf2 = nanmean(tLists{2}.f2,1);
        
        %normalize rivalry transitions using sim amplitudes
        l2hRivf1 = (l2hRivf1 - min(l2hRivf1)) / f1simAmp;
        h2lRivf1 = (h2lRivf1 - min(h2lRivf1)) / f1simAmp;
        l2hRivf2 = (l2hRivf2 - min(l2hRivf2)) / f2simAmp;
        h2lRivf2 = (h2lRivf2 - min(h2lRivf2)) / f2simAmp;
        
        % normalize rivalry transitions using sim transitions
        if size(tLists{1}.f1,1) > minInstances
            groupTLists{7}.f1 = [groupTLists{7}.f1; l2hRivf1];
            groupTLists{7}.f2 = [groupTLists{7}.f2; l2hRivf2];
        end
        if size(tLists{2}.f1,1) > minInstances
            groupTLists{8}.f1 = [groupTLists{8}.f1; h2lRivf1];
            groupTLists{8}.f2 = [groupTLists{8}.f2; h2lRivf2];
        end
        
       
        %% Collect transitions
        for tType = 1:length(tLists)
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
            end
        end
        
        %% Collect dartRival peaks normalized to dartSim amplitudes
        
        % Find average sim amplitude
        
        %collect sim maxes and mins       
        for pType = 3:4
            if pType == 3 % 3: sim low dom
                f1lowMaxes = max(pLists{pType}.f1, [],2);
                f1lowMins = min(pLists{pType}.f1, [], 2);
                f2lowMaxes = max(pLists{pType}.f2, [],2);
                f2lowMins = min(pLists{pType}.f2, [], 2);
            elseif pType == 4 % 4: sim high dom
                f1highMaxes = max(pLists{pType}.f1, [],2);
                f1highMins = min(pLists{pType}.f1, [], 2);
                f2highMaxes = max(pLists{pType}.f2, [],2);
                f2highMins = min(pLists{pType}.f2, [], 2);
            end
        end
        
        % calculate average peak to trough sim amplitude
        f1lsimAmp = nanmean(f1lowMaxes) - nanmean(f1lowMins);
        f1hsimAmp = nanmean(f1highMaxes) - nanmean(f1highMins);
        f2lsimAmp = nanmean(f2lowMaxes) - nanmean(f2lowMins);
        f2hsimAmp = nanmean(f2highMaxes) - nanmean(f2highMins);
        
        % mean of participant's rivalry peaks
        ldomRivf1 = nanmean(pLists{1}.f1,1);
        hdomRivf1 = nanmean(pLists{2}.f1,1);
        ldomRivf2 = nanmean(pLists{1}.f2,1);
        hdomRivf2 = nanmean(pLists{2}.f2,1);
        
        %normalize rivalry peaks using sim amplitudes
        ldomRivf1 = (ldomRivf1 - min(ldomRivf1)) / f1lsimAmp;
        hdomRivf1 = (hdomRivf1 - min(hdomRivf1)) / f1hsimAmp;
        ldomRivf2 = (ldomRivf2 - min(ldomRivf2)) / f2lsimAmp;
        hdomRivf2 = (hdomRivf2 - min(hdomRivf2)) / f2hsimAmp;
        
        % collect peaks
        if size(pLists{1}.f1,1) > minInstances
            groupPLists{7}.f1 = [groupPLists{7}.f1; ldomRivf1];
            groupPLists{7}.f2 = [groupPLists{7}.f2; ldomRivf2];
        end
        if size(pLists{2}.f1,1) > minInstances
            groupPLists{8}.f1 = [groupPLists{8}.f1; hdomRivf1];
            groupPLists{8}.f2 = [groupPLists{8}.f2; hdomRivf2];
        end
        
        %normAmp = mean([max(l2hRivf1) max(h2lRivf1) max(l2hRivf2) max(h2lRivf2)]);
        %save([indicesDir 'normAmps/' parName], 'normAmp')
        
        %% Collect peaks
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
            end
        end
    end
    
    %% Plotting
    
    % transitions
    %for tType = 1:length(groupTLists)
    for tType = 7:8
        
        % must have a minimum number of transitions for each run type
        if size(groupTLists{tType}.f1, 1) < minPars
            disp([analGroupIDs{iGroup} 'does not have enough ' tLabels(tType) ' participants'])
            continue;
        end
    
        meanF1trace = nanmean(groupTLists{tType}.f1,1);
        errorF1trace = ste(groupTLists{tType}.f1);
        
        meanF2trace = nanmean(groupTLists{tType}.f2,1);
        errorF2trace = ste(groupTLists{tType}.f2);
        
        figure
        title([analGroupIDs{iGroup} ' ' tLabels{tType} ': ' num2str(size(groupTLists{tType}.f1, 1)) ' participants averaged'])
        hold on
        mseb([(-transHalf):(1/512):transHalf],[meanF1trace; meanF2trace],[errorF1trace; errorF2trace],[],1);
        if tType == 7 || tType == 8
            ylim([0 1])
        else
            ylim([1 2])
        end
        vline(0,'k','transition')
        xlabel('Time from button press (s)');
        ylabel('Amplitude');
        legend('Low Frequency', 'High Frequency');
        
        if allTrans
            figName = [groupPlotDir analGroupIDs{iGroup} '_' tFileNames{tType} '_' num2str(size(groupTLists{tType}.f1,1)) 'pars_ALL.jpg'];
        else
            figName = [groupPlotDir analGroupIDs{iGroup} '_' tFileNames{tType} '_' num2str(size(groupTLists{tType}.f1,1)) 'pars.jpg'];
        end
        saveas(gcf,figName,'jpg')
    end
    
    % peaks
     %for pType = 1:length(groupPLists)
     for pType = 7:8
        
        % must have a minimum number of transitions for each run type
        if size(groupPLists{pType}.f1, 1) < minInstances
            disp([analGroupIDs{iGroup} 'does not have enough ' pLabels(pType) ' instances'])
            continue;
        end
    
        meanF1trace = nanmean(groupPLists{pType}.f1,1);
        errorF1trace = ste(groupPLists{pType}.f1);
        
        meanF2trace = nanmean(groupPLists{pType}.f2,1);
        errorF2trace = ste(groupPLists{pType}.f2);
        
        figure
        title([analGroupIDs{iGroup} ' ' pLabels{pType} ': ' num2str(size(groupPLists{pType}.f1, 1)) ' participants averaged'])
        hold on
        mseb([(-peakHalf):(1/512):peakHalf],[meanF1trace; meanF2trace],[errorF1trace; errorF2trace],[],1);
        if pType == 7 || pType == 8
            ylim([0 1])
        else
            ylim([1 1.9])
        end
        vline(0,'k','peak')
        xlabel('Time from peak (s)');
        ylabel('Amplitude');
        legend('Low Frequency', 'High Frequency');
        
        if altPeak
            figName = [groupPlotDir analGroupIDs{iGroup} '_' pFileNames{pType} 'Peaks_' num2str(size(groupPLists{pType}.f1,1)) 'pars_ALT.jpg'];
        else
        figName = [groupPlotDir analGroupIDs{iGroup} '_' pFileNames{pType} 'Peaks_' num2str(size(groupPLists{pType}.f1,1)) 'pars.jpg'];
        end
        saveas(gcf,figName,'jpg')
    end
end
end


