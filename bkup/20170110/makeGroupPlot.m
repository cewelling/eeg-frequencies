function makeGroupPlot()
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% load parameters
groupAnalysisParams

numFormat = '%02d'; % for participant number

tLabels = {'dartboard rivalry low to high' 'dartboard rivalry high to low' 'sim low to high' 'sim high to low' ...
    'marz rivalry low to high' 'marz rivalry high to low'};
tFileNames = {'dartRival_l2h' 'dartRival_h2l' 'dartSim_l2h' 'dartSim_h2l' 'marzRival_l2h' 'marzRival_h2l'};

pLabels = {'dartboard rivalry low dom' 'dartboard rivalry high dom' 'sim low dom' 'sim high dom' ...
    'marz rivalry low dom' 'marz rivalry high dom'};
pFileNames = {'dartRival_lowDom' 'dartRival_highDom' 'simRival_lowDom' 'simRival_highDom' ...
    'marzRival_lowDom' 'marzRival_highDom'};

for iGroup = 1:length(groupParNums)
    
    % allocate space to store transitions
    l2hRival.f1 = []; l2hRival.f2 = [];
    h2lRival.f1 = []; h2lRival.f2 = [];
    l2hSim.f1 = []; l2hSim.f2 = [];
    h2lSim.f1 = []; h2lSim.f2 = [];
    l2hMarzRival.f1 = []; l2hMarzRival.f2 = [];
    h2lMarzRival.f1 = []; h2lMarzRival.f2 = [];
    groupTLists = {l2hRival h2lRival l2hSim h2lSim l2hMarzRival h2lMarzRival};
    
    % allocate space to store peaks
    slowPkRival.f1 = []; slowPkRival.f2 = [];
    fastPkRival.f1 = []; fastPkRival.f2 = [];
    slowPkSim.f1 = []; slowPkSim.f2 = [];
    fastPkSim.f1 = []; fastPkSim.f2 = [];
    slowPkMarzRival.f1 = []; slowPkMarzRival.f2 = [];
    fastPkMarzRival.f1 = []; fastPkMarzRival.f2 = [];
    groupPLists = {slowPkRival fastPkRival slowPkSim fastPkSim slowPkMarzRival fastPkMarzRival};
    
    for parNum = groupParNums{iGroup};
        parName = [groupCodes{iGroup} num2str(parNum, numFormat)];
        
        % load participant's transitions (if participant exists)
        if allTrans
            if exist(['transitions2/' parName '.mat'], 'file')
                load(['transitions2/' parName '.mat'])
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
        if exist(['peaks/' parName '.mat'], 'file')
            load(['peaks/' parName '.mat'])
        else
            continue;
        end
        
        % Collect transitions
        for tType = 1:length(tLists)
            
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
        
        % Collect peaks
        for pType = 1:length(pLists)
            
            % mean of participant's transitions
            f1mean = nanmean(pLists{pType}.f1,1);
            f2mean = nanmean(pLists{pType}.f2,1);
            if isnan(f1mean)
                f1mean = [];
            end
            if isnan(f2mean)
                f2mean = [];
            end
            
            % store transitions
            groupPLists{pType}.f1 = [groupPLists{pType}.f1; f1mean];
            groupPLists{pType}.f2 = [groupPLists{pType}.f2; f2mean];
        end
    end
    
    % Plotting
    
    % transitions
    for tType = 1:length(groupTLists)
        
        % must have a minimum number of transitions for each run type
        if size(groupTLists{tType}.f1, 1) < minInstances
            disp([analGroupIDs{iGroup} 'does not have enough ' tLabels(tType) ' instances'])
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
        ylim([-0.6 0.8])
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
     for pType = 1:length(groupPLists)
        
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
        ylim([-0.4 1.2])
        vline(0,'k','peak')
        xlabel('Time from peak (s)');
        ylabel('Amplitude');
        legend('Low Frequency', 'High Frequency');
        
        
        figName = [groupPlotDir analGroupIDs{iGroup} '_' pFileNames{pType} 'Peaks_' num2str(size(groupPLists{pType}.f1,1)) 'pars.jpg'];
        saveas(gcf,figName,'jpg')
    end
end
end


