function makeGroupPlot()
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% load parameters
groupAnalysisParams

numFormat = '%02d'; % for participant number

tLabels = {'dartboard rivalry low to high' 'dartboard rivalry high to low' 'sim low to high' 'sim high to low' ...
    'marz rivalry low to high' 'marz rivalry high to low'};

for iGroup = 1:length(groupParNums)
    
    % allocate space to store transitions
    l2hRival.f1 = []; l2hRival.f2 = [];
    h2lRival.f1 = []; h2lRival.f2 = [];
    l2hSim.f1 = []; l2hSim.f2 = [];
    h2lSim.f1 = []; h2lSim.f2 = [];
    l2hMarzRival.f1 = []; l2hMarzRival.f2 = [];
    h2lMarzRival.f1 = []; h2lMarzRival.f2 = [];
    groupTLists = {l2hRival h2lRival l2hSim h2lSim l2hMarzRival h2lMarzRival};
    
    for parNum = groupParNums{iGroup};
        parName = [groupCodes{iGroup} num2str(parNum, numFormat)];
        
        % load participant's transitions (if participant exists)
        if exist(['transitions/' parName '.mat'], 'file')
            load(['transitions/' parName '.mat'])
        else
            continue;
        end
        
        % Collect transitions
        for tType = 1:length(tLists)
            groupTLists{tType}.f1 = [groupTLists{tType}.f1; tLists{tType}.f1];
            groupTLists{tType}.f2 = [groupTLists{tType}.f2; tLists{tType}.f2];
        end
    end
    
    % Plotting
    
    for tType = 1:length(groupTLists)
        
        % must have a minimum number of transitions for each run type
        if size(groupTLists{tType}.f1, 1) < minInstances
            disp([analGroupIDs{iGroup} 'does not have enough ' tLabels(tType) ' instances'])
            continue;
        end
    
        meanF1trace = mean(groupTLists{tType}.f1,1);
        errorF1trace = ste(groupTLists{tType}.f1);
        
        meanF2trace = mean(groupTLists{tType}.f2,1);
        errorF2trace = ste(groupTLists{tType}.f2);
        
        figure
        title([analGroupIDs{iGroup} ' ' tLabels{tType} ': ' num2str(size(groupTLists{tType}.f1, 1)) ' transitions averaged'])
        hold on
        mseb([(-transHalf):(1/512):transHalf],[meanF1trace; meanF2trace],[errorF1trace; errorF2trace],[],1);
        vline(0,'k','transition')
        xlabel('Time from button press (s)');
        ylabel('Amplitude');
        legend('Low Frequency', 'High Frequency');
    end
end
end


