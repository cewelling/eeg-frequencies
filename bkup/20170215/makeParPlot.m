function [ tLists, tLabels, pLists, pLabels] = makeParPlot( parName, date)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%import parameters
analysisParams

% only take runs above min SNR
runList = {};
excluded = {};
for iRunType = 1:length(runTypes)
    cRunType = runTypes{iRunType};
    
    % load or (if not yet set) set common electrodes for each type of run based on SNR
    if strfind(cRunType, 'dart')
        if exist(['highSNRelecs/' parName '_darts_' num2str(commonElecs) 'of' num2str(numElecs) 'elecs.mat'], 'file')
            load(['highSNRelecs/' parName '_darts_' num2str(commonElecs) 'of' num2str(numElecs) 'elecs']);
        else
            elecs = setElecs(parName);
        end
    else
        if exist(['highSNRelecs/' parName '_marz_' num2str(commonElecs) 'of' num2str(numElecs) 'elecs.mat'], 'file')
            load(['highSNRelecs/' parName '_marz_' num2str(commonElecs) 'of' num2str(numElecs) 'elecs']);
        else
            elecs = setElecs(parName);
        end
    end
    
    for runIndex = runIndices
        runName = [cRunType num2str(runIndex)];
        
        % Get the appropriate EEG file for this run
        EEGfile = [eegDir parName '_' runName '_' date '.bdf'];
        
        % skip over participants / runs that don't exist
        if ~exist(EEGfile, 'file')
            continue;
        end
        
        % Load SNRs of high SNR electrodes for this run
        load(['highSNRelecs/' parName '_' runName '_'  num2str(numElecs) 'elecs.mat']);
        if nanmean(maxSNRs(2,ismember(maxSNRs(1,:), elecs))) > minSNR
            runList = [runList runName];
        else
            excluded = [excluded runName];
        end
    end
end

% save list of runs excluded due to low SNR
save(['runsExcluded/' parName], 'excluded');

% transition case codes and peak types
l2h = 1; %low to high
h2l = 2;
l2l = 3;
h2h = 4;

slow = 1; % peaks in slow freq series
fast = 2;

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

% allocate space to store peaks
slowPkRival.f1 = []; slowPkRival.f2 = [];
fastPkRival.f1 = []; fastPkRival.f2 = [];
slowPkSim.f1 = []; slowPkSim.f2 = [];
fastPkSim.f1 = []; fastPkSim.f2 = [];
slowPkMarzRival.f1 = []; slowPkMarzRival.f2 = [];
fastPkMarzRival.f1 = []; fastPkMarzRival.f2 = [];

% collect transitions and peaks, keep rival types and sim separate
for i = 1:length(runList)
    runName = runList{i};
    
    % Get the appropriate EEG file for this run
    EEGfile = [eegDir parName '_' runName '_' date '.bdf'];
    
    % load transitions for this run
     if ~exist(['tListsRLS/' parName '_' runName '_' normType '.mat'], 'file')
        transitions = aggTransitions(parName, runName, date, EEGfile);
     else
         load(['tListsRLS/' parName '_' runName '_' normType '.mat'])
     end
    
    % load peaks for this run
    if altPeak &&  ~exist(['altpListsRLS/' parName '_' runName '_' normType '.mat'], 'file') || ...
            ~altPeak && ~exist(['pListsRLS/' parName '_' runName '_' normType '.mat'], 'file')
        peaks = aggPeaks(parName, runName, date, EEGfile);
    elseif altPeak &&  exist(['altpListsRLS/' parName '_' runName '_' normType '.mat'], 'file')
        load(['altpListsRLS/' parName '_' runName '_' normType '.mat'])
    else
        load(['pListsRLS/' parName '_' runName '_' normType '.mat'])
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
            if criteria(i,1) >= domMin && criteria(i,3) <= gapMax
                if ~noMixed && (criteria(i,2) >= mixedMin || isnan(criteria(i,2)))
                    tKeep{cCode} = [tKeep{cCode} i];
                elseif noMixed && isnan(criteria(i,2))
                    tKeep{cCode} = [tKeep{cCode} i];
                end
            end
        end
    end
    
    % determine which peaks fit criteria
    pKeep = [];
    for domFreq = 1:2 % low, high
        criteria = peaks(domFreq).durs;
        pKeep{domFreq} = [];
        if isempty(criteria)
            continue;
        end
        for i = 1:size(criteria, 1)
            if criteria(i) > minPwidth && criteria(i) < maxPwidth
                pKeep{domFreq} = [pKeep{domFreq} i];
            end
        end
    end
    
    if strfind(runName, 'dartRival')
        % dartboard rival transitions
        l2hRival.f1 = [l2hRival.f1; transitions(l2h).f1(tKeep{1}, :)];
        l2hRival.f2 = [l2hRival.f2; transitions(l2h).f2(tKeep{1}, :)];
        h2lRival.f1 = [h2lRival.f1; transitions(h2l).f1(tKeep{2}, :)];
        h2lRival.f2 = [h2lRival.f2; transitions(h2l).f2(tKeep{2}, :)];
        l2lRival.f1 = [l2lRival.f1; transitions(l2l).f1(tKeep{3}, :)];
        l2lRival.f2 = [l2lRival.f2; transitions(l2l).f2(tKeep{3}, :)];
        h2hRival.f1 = [h2hRival.f1; transitions(h2h).f1(tKeep{4}, :)];
        h2hRival.f2 = [h2hRival.f2; transitions(h2h).f2(tKeep{4}, :)];
        
        % dartboard rival peaks
        slowPkRival.f1 = [slowPkRival.f1; peaks(slow).f1(pKeep{1}, :)];
        slowPkRival.f2 = [slowPkRival.f2; peaks(slow).f2(pKeep{1}, :)];
        fastPkRival.f1 = [fastPkRival.f1; peaks(fast).f1(pKeep{2}, :)];
        fastPkRival.f2 = [fastPkRival.f2; peaks(fast).f2(pKeep{2}, :)];
        
    elseif strfind(runName, 'dartSim')
        % dartboard sim transitions
        l2hSim.f1 = [l2hSim.f1; transitions(l2h).f1(tKeep{1}, :)];
        l2hSim.f2 = [l2hSim.f2; transitions(l2h).f2(tKeep{1}, :)];
        h2lSim.f1 = [h2lSim.f1; transitions(h2l).f1(tKeep{2}, :)];
        h2lSim.f2 = [h2lSim.f2; transitions(h2l).f2(tKeep{2}, :)];
        
        % dartboard sim peaks
        slowPkSim.f1 = [slowPkSim.f1; peaks(slow).f1(pKeep{1}, :)];
        slowPkSim.f2 = [slowPkSim.f2; peaks(slow).f2(pKeep{1}, :)];
        fastPkSim.f1 = [fastPkSim.f1; peaks(fast).f1(pKeep{2}, :)];
        fastPkSim.f2 = [fastPkSim.f2; peaks(fast).f2(pKeep{2}, :)];
        
    elseif strfind(runName, 'marzRival')
        % marz rival transitions
        l2hMarzRival.f1 = [l2hRival.f1; transitions(l2h).f1(tKeep{1}, :)];
        l2hMarzRival.f2 = [l2hRival.f2; transitions(l2h).f2(tKeep{1}, :)];
        h2lMarzRival.f1 = [h2lRival.f1; transitions(h2l).f1(tKeep{2}, :)];
        h2lMarzRival.f2 = [h2lRival.f2; transitions(h2l).f2(tKeep{2}, :)];
        l2lMarzRival.f1 = [l2lRival.f1; transitions(l2l).f1(tKeep{3}, :)];
        l2lMarzRival.f2 = [l2lRival.f2; transitions(l2l).f2(tKeep{3}, :)];
        h2hMarzRival.f1 = [h2hRival.f1; transitions(h2h).f1(tKeep{4}, :)];
        h2hMarzRival.f2 = [h2hRival.f2; transitions(h2h).f2(tKeep{4}, :)];
        
        % marz rival peaks
        slowPkMarzRival.f1 = [slowPkMarzRival.f1; peaks(slow).f1(pKeep{1}, :)];
        slowPkMarzRival.f2 = [slowPkMarzRival.f2; peaks(slow).f2(pKeep{1}, :)];
        fastPkMarzRival.f1 = [fastPkMarzRival.f1; peaks(fast).f1(pKeep{2}, :)];
        fastPkMarzRival.f2 = [fastPkMarzRival.f2; peaks(fast).f2(pKeep{2}, :)];
    end
end

tLists = {l2hRival h2lRival l2lRival h2hRival l2hSim h2lSim l2hMarzRival h2lMarzRival l2lMarzRival h2hMarzRival};
tLabels = {'dartboard rivalry low to high' 'dartboard rivalry high to low' ...
    'dartboard rivalry low to low' 'dartboard rivalry high to high' ...
    'sim low to high' 'sim high to low' 'marz rivalry low to high' 'marz rivalry high to low' ...
    'marz rivalry low to low' 'marz rivalry high to high'};

save (['transitions/' parName '_' normType], 'tLists');

pLists = {slowPkRival fastPkRival slowPkSim fastPkSim slowPkMarzRival fastPkMarzRival};
pLabels = {'dartboard rivalry low dom' 'dartboard rivalry high dom' 'sim low dom' 'sim high dom' ...
    'marz rivalry low dom' 'marz rivalry high dom'};
pFileNames = {'dartRival_lowDom' 'dartRival_highDom' 'simRival_lowDom' 'simRival_highDom' ...
    'marzRival_lowDom' 'marzRival_highDom'};

save(['peaks/' parName '_' normType], 'pLists');

if strcmp(parPlotOrNot, 'yes')
    %plot averaged transitions
    for tType = 1:length(tLists)
    %for tType = 5:6    
        % must have a minimum number of transitions for each run type
        if size(tLists{tType}.f1, 1) < minInstances
            disp([parName 'does not have enough ' tLabels(tType) ' instances'])
            continue;
        end
        
        meanF1trace = nanmean(tLists{tType}.f1,1);
        errorF1trace = ste(tLists{tType}.f1);
        
        meanF2trace = nanmean(tLists{tType}.f2,1);
        errorF2trace = ste(tLists{tType}.f2);
        
        figure
        title([parName ' ' tLabels{tType} ': ' num2str(size(tLists{tType}.f1, 1)) ' transitions averaged'])
        hold on
        mseb([(-transHalf):(1/512):transHalf],[meanF1trace; meanF2trace],[errorF1trace; errorF2trace],[],1);
        vline(0,'k','transition')
        xlabel('Time from button press (s)');
        ylabel('Amplitude');
        legend('Low Frequency', 'High Frequency');
    end
    
%     % plot averaged peaks
    %for pType = 1:length(pLists)
%     for pType = 1:2
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
%         title([parName ' ' pLabels{pType} ': ' num2str(size(pLists{pType}.f1, 1)) ' peaks averaged'])
%         hold on
%         mseb([(-peakHalf):(1/512):peakHalf],[meanF1trace; meanF2trace],[errorF1trace; errorF2trace],[],1);
%         vline(0,'k','peak')
%         xlabel('Time from peak (s)');
%         ylabel('Amplitude');
%         legend('Low Frequency', 'High Frequency');
%         
%         %disp(['Runs included: ' runList])
%         
%         % end
%         
%         %lowDomAmp
%         
%         % Calculate and save modulation index
%         lowDomMod = nanmax(nanmean(pLists{1}.f1, 1)) / nanmax(nanmean(pLists{3}.f1, 1));
%         highDomMod = nanmax(nanmean(pLists{2}.f2, 1)) / nanmax(nanmean(pLists{4}.f2, 1));
%         save([indicesDir 'modIndices/' parName '_lowDom'], 'lowDomMod');
%         save([indicesDir 'modIndices/' parName '_highDom'], 'highDomMod');
%         
%     end
end

end


