function [ tLists tLabels pLists pLabels] = makeParPlot( parName, date, EEGfiles, iPar )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%import parameters
analysisParams  

% set common electrodes for each type of run based on SNR
setElecs(parName)

% using those electrodes, only take runs above min SNR
runList = {};
for iRunType = 1:length(runTypes)
    cRunType = runTypes{iRunType};
    for runIndex = runIndices
        runName = [cRunType num2str(runIndex)];
        
        % Get the appropriate EEG file for this run
        EEGfile = [eegDir parName '_' runName '_' date '.bdf'];
        
        % skip over participants / runs that don't exist
        if ~exist(EEGfile, 'file')
            continue;
        end
        
        % Load
        if strfind(cRunType, 'dart')
            load(['highSNRelecs/' parName '_darts_' num2str(commonElecs) 'of' num2str(numElecs) 'elecs']);
        else
            load(['highSNRelecs/' parName '_marz_' num2str(commonElecs) 'of' num2str(numElecs) 'elecs']);
        end            
        
        % Load SNRs of high SNR electrodes for this run
        load(['highSNRelecs/' parName '_' runName '_'  num2str(numElecs) 'elecs.mat']);
        if nanmean(maxSNRs(2,ismember(maxSNRs(1,:), elecs))) > minSNR
            runList = [runList runName];
        end
    end
end

% transition case codes and peak types
l2h = 1; %low to high
h2l = 2;
slow = 1; % peaks in slow freq series
fast = 2;

% allocate space to store transitions
l2hRival.f1 = []; l2hRival.f2 = [];
h2lRival.f1 = []; h2lRival.f2 = [];
l2hSim.f1 = []; l2hSim.f2 = [];
h2lSim.f1 = []; h2lSim.f2 = [];
l2hMarzRival.f1 = []; l2hMarzRival.f2 = [];
h2lMarzRival.f1 = []; h2lMarzRival.f2 = [];

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
    transitions = aggTransitions(parName, runName, date, EEGfile);
    peaks = aggPeaks(parName, runName, date, EEGfile);
    
    if strfind(runName, 'dartRival')
        % dartboard rival transitions
        l2hRival.f1 = [l2hRival.f1; transitions(l2h).f1];
        l2hRival.f2 = [l2hRival.f2; transitions(l2h).f2];
        h2lRival.f1 = [h2lRival.f1; transitions(h2l).f1];
        h2lRival.f2 = [h2lRival.f2; transitions(h2l).f2];
        
        % dartboard rival peaks
        slowPkRival.f1 = [slowPkRival.f1; peaks(slow).f1];
        slowPkRival.f2 = [slowPkRival.f2; peaks(slow).f2];
        fastPkRival.f1 = [fastPkRival.f1; peaks(fast).f1];
        fastPkRival.f2 = [fastPkRival.f2; peaks(fast).f2];
        
    elseif strfind(runName, 'dartSim')
        % dartboard sim transitions
        l2hSim.f1 = [l2hSim.f1; transitions(l2h).f1];
        l2hSim.f2 = [l2hSim.f2; transitions(l2h).f2];
        h2lSim.f1 = [h2lSim.f1; transitions(h2l).f1];
        h2lSim.f2 = [h2lSim.f2; transitions(h2l).f2];
        
        % dartboard sim peaks
        slowPkSim.f1 = [slowPkSim.f1; peaks(slow).f1];
        slowPkSim.f2 = [slowPkSim.f2; peaks(slow).f2];
        fastPkSim.f1 = [fastPkSim.f1; peaks(fast).f1];
        fastPkSim.f2 = [fastPkSim.f2; peaks(fast).f2];

    elseif strfind(runName, 'marzRival')
        % marz rival transitions
        l2hMarzRival.f1 = [l2hRival.f1; transitions(l2h).f1];
        l2hMarzRival.f2 = [l2hRival.f2; transitions(l2h).f2];
        h2lMarzRival.f1 = [h2lRival.f1; transitions(h2l).f1];
        h2lMarzRival.f2 = [h2lRival.f2; transitions(h2l).f2];
        
        % marz rival peaks
        slowPkMarzRival.f1 = [slowPkMarzRival.f1; peaks(slow).f1];
        slowPkMarzRival.f2 = [slowPkMarzRival.f2; peaks(slow).f2];
        fastPkMarzRival.f1 = [fastPkMarzRival.f1; peaks(fast).f1];
        fastPkMarzRival.f2 = [fastPkMarzRival.f2; peaks(fast).f2];
    end
end

tLists = {l2hRival h2lRival l2hSim h2lSim l2hMarzRival h2lMarzRival};
tLabels = {'dartboard rivalry low to high' 'dartboard rivalry high to low' 'sim low to high' 'sim high to low' ...
    'marz rivalry low to high' 'marz rivalry high to low'};
if allTrans
    if strcmp(normType, 'none')
        save(['transNoNorm2/' parName], 'tLists');
    else
        save(['transitions2/' parName], 'tLists');
    end
elseif strcmp(normType, 'none')
    save(['transNoNorm/' parName], 'tLists');
else
    save (['transitions/' parName], 'tLists');
end

pLists = {slowPkRival fastPkRival slowPkSim fastPkSim slowPkMarzRival fastPkMarzRival};
pLabels = {'dartboard rivalry low dom' 'dartboard rivalry high dom' 'sim low dom' 'sim high dom' ...
    'marz rivalry low dom' 'marz rivalry high dom'};
pFileNames = {'dartRival_lowDom' 'dartRival_highDom' 'simRival_lowDom' 'simRival_highDom' ...
    'marzRival_lowDom' 'marzRival_highDom'};

if altPeak
    if strcmp(normType, 'none')
        save(['peaksNoNorm2/' parName], 'pLists');
    else
        save(['peaks2/' parName], 'pLists');
    end
elseif strcmp(normType, 'none')
    save(['peaksNoNorm/' parName], 'pLists');
else
    save (['peaks/' parName], 'pLists');
end

%plot averaged transitions
load(['transNoNorm2/' parName '.mat'])
for tType = 1:length(tLists)

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

% plot averaged peaks
load(['peaks/' parName '.mat'])
for pType = 1:length(pLists)
    
    % must have a minimum number of peaks for each run type
    if size(pLists{pType}.f1, 1) < minInstances
        disp([parName 'does not have enough ' pLabels(pType) ' peaks'])
        continue;
    end
    
    meanF1trace = nanmean(pLists{pType}.f1,1);
    errorF1trace = ste(pLists{pType}.f1);
    
    meanF2trace = nanmean(pLists{pType}.f2,1);
    errorF2trace = ste(pLists{pType}.f2);
    
    figure
    title([parName ' ' pLabels{pType} ': ' num2str(size(pLists{pType}.f1, 1)) ' peaks averaged'])
    hold on
    mseb([(-peakHalf):(1/512):peakHalf],[meanF1trace; meanF2trace],[errorF1trace; errorF2trace],[],1);
    vline(0,'k','peak')
    xlabel('Time from peak (s)');
    ylabel('Amplitude');
    legend('Low Frequency', 'High Frequency');

%disp(['Runs included: ' runList])

% end
    
%lowDomAmp 

% Calculate and save modulation index
lowDomMod = nanmax(nanmean(pLists{1}.f1, 1)) / nanmax(nanmean(pLists{3}.f1, 1));
highDomMod = nanmax(nanmean(pLists{2}.f2, 1)) / nanmax(nanmean(pLists{4}.f2, 1));
save([indicesDir 'modIndices/' parName '_lowDom'], 'lowDomMod');
save([indicesDir 'modIndices/' parName '_highDom'], 'highDomMod');

end


