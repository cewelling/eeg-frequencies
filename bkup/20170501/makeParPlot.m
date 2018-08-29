function [ tLists, tLabels, pLists, pLabels] = makeParPlot( parName, date)
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

%% collect all transitions and peaks for each participant

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
        
        % load transitions for this run
        if exist(['transitions/tListsRLS/' parName '_' runName '_' normType '.mat'], 'file') && isequal(analFreqs,stimFreqs)
            load(['transitions/tListsRLS/' parName '_' runName '_' normType '.mat'])
        elseif exist(['transitions/tListsRLS/' parName '_' runName '_' normType '_harmonics.mat'], 'file') && isequal(analFreqs,harFreqs)
            load(['transitions/tListsRLS/' parName '_' runName '_' normType '_harmonics.mat'])
        else
            transitions = aggTransitions(parName, runName, date, EEGfile);
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
        
        % determine which transitions fit criteria
        tKeep = [];
        for cCode = 1:4 % l2h, h2l, l2l, h2h
            criteria = transitions(cCode).crit;
            tKeep{cCode} = [];
            if isempty(criteria)
                continue;
            end
            for i = 1:size(criteria, 1)
                if criteria(i,6) >= minSNR % exclude transitions by SNR
                    if criteria(i,1) >= dom1Min && criteria(i,1) <= dom1Max && criteria(i,4) <= gapMax
                        if criteria(i,2) >= dom2Min && criteria(i,2) < dom2Max
                            if ~noMixed && ((criteria(i,3) >= mixedMin && criteria(i,3) <= mixedMax))
                                tKeep{cCode} = [tKeep{cCode} i];
                            elseif noMixed && criteria(i,3) == 0
                                tKeep{cCode} = [tKeep{cCode} i];
                            end
                        end
                    end
                end
            end
        end
        
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
                    if peakSNRs(i) > minSNR
                        pKeep{domFreq} = [pKeep{domFreq} i];
                    end
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
            l2hSim.f1 = [l2hSim.f1; transitions(l2h).f1(:, :)]; % no criteria for simulations
            l2hSim.f2 = [l2hSim.f2; transitions(l2h).f2(:, :)];
            h2lSim.f1 = [h2lSim.f1; transitions(h2l).f1(:, :)];
            h2lSim.f2 = [h2lSim.f2; transitions(h2l).f2(:, :)];
            
            % dartboard sim peaks
            slowPkSim.f1 = [slowPkSim.f1; peaks(slow).f1(:, :)];
            slowPkSim.f2 = [slowPkSim.f2; peaks(slow).f2(:, :)];
            fastPkSim.f1 = [fastPkSim.f1; peaks(fast).f1(:, :)];
            fastPkSim.f2 = [fastPkSim.f2; peaks(fast).f2(:, :)];
            
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
end

%% Save collected transitions and peaks

tLists = {l2hRival h2lRival l2lRival h2hRival l2hSim h2lSim l2hMarzRival h2lMarzRival l2lMarzRival h2hMarzRival};
tLabels = {'dartboard rivalry low to high' 'dartboard rivalry high to low' ...
    'dartboard rivalry low to low' 'dartboard rivalry high to high' ...
    'sim low to high' 'sim high to low' 'marz rivalry low to high' 'marz rivalry high to low' ...
    'marz rivalry low to low' 'marz rivalry high to high'};

save (['transitions/finalTransitions' parName '_' normType], 'tLists');

pLists = {slowPkRival fastPkRival slowPkSim fastPkSim slowPkMarzRival fastPkMarzRival};
pLabels = {'dartboard rivalry low dom' 'dartboard rivalry high dom' 'sim low dom' 'sim high dom' ...
    'marz rivalry low dom' 'marz rivalry high dom'};
pFileNames = {'dartRival_lowDom' 'dartRival_highDom' 'simRival_lowDom' 'simRival_highDom' ...
    'marzRival_lowDom' 'marzRival_highDom'};

save(['peaks/finalPeaks' parName '_' normType], 'pLists');

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
    
    % must have a minimum number of transitions
    if size([tLists{1}.f2;tLists{2}.f1]) >= minInstances
        
        % plot dartRival transitions: l2h and h2l averaged
        s2dTrace = nanmean([tLists{1}.f2; tLists{2}.f1],1);
        s2dError = ste([tLists{1}.f2; tLists{2}.f1]);
        d2sTrace = nanmean([tLists{1}.f1; tLists{2}.f2],1);
        d2sError = ste([tLists{1}.f1; tLists{2}.f2]);
        
        figure
        title([parName ' dartRival: ' num2str(size([tLists{1}.f1; tLists{2}.f1], 1)) ' transitions averaged'])
        hold on
        mseb([(-transHalf):(1/sampRate):transHalf],[s2dTrace; d2sTrace],[s2dError; d2sError],[],1);
        ylim([-1.5 1.5])
        vline(0,'k','transition')
        xlabel('Time from button press (s)');
        ylabel('Amplitude');
        legend('Frequency 1', 'Frequency 2');       
    else
        disp([parName 'does not have enough transitions'])
    end
    
    %% plot averaged peaks
    
    for pType = 1:length(pLists) % 1:2
        
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
    end
    %         % Calculate and save modulation index
    %         lowDomMod = nanmax(nanmean(pLists{1}.f1, 1)) / nanmax(nanmean(pLists{3}.f1, 1));
    %         highDomMod = nanmax(nanmean(pLists{2}.f2, 1)) / nanmax(nanmean(pLists{4}.f2, 1));
    %         save([indicesDir 'modIndices/' parName '_lowDom'], 'lowDomMod');
    %         save([indicesDir 'modIndices/' parName '_highDom'], 'highDomMod');
end

end


