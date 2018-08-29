function [ modIndex ] = getModIndex(parName, date, EEGfiles, iPar)
%UNTITLED18 Summary of this function goes here
%   Detailed explanation goes here

% import parameters
analysisParams

% only take runs above min SNR
runList = [];
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
        
        SNRs = runFFT(EEGfile, runName, parName, iPar);
        if length(SNRs) > 2
            error('Choose just one noise window, electrode set for %s run %d SNR', parName, runName)
        end
        if mean([SNRs.value]) > 3
            runList = [runList runName];
        end
    end
end

% case codes
l2h = 1; %low to high
h2l = 2;

% allocate space to store transitions
l2hRival.f1 = []; l2hRival.f2 = [];
h2lRival.f1 = []; h2lRival.f2 = [];
l2hSim.f1 = []; l2hSim.f2 = [];
h2lSim.f1 = []; h2lSim.f2 = [];

% collect transitions, keep rival and sim separate
for i = 1:length(runList)
    runName = runList{i};
    % Get the appropriate EEG file for this run
    EEGfile = [eegDir parName '_' runName '_' date '.bdf'];
    transitions = aggTransitions(parName, runName, date, EEGfile);
    if strfind(runName, 'dartRival')
        l2hRival.f1 = [l2hRival.f1; transitions(l2h).f1];
        l2hRival.f2 = [l2hRival.f2; transitions(l2h).f2];
        h2lRival.f1 = [h2lRival.f1; transitions(h2l).f1];
        h2lRival.f2 = [h2lRival.f2; transitions(h2l).f2];
    elseif strfind(runName, 'dartSim')
        l2hSim.f1 = [l2hSim.f1; transitions(l2h).f1];
        l2hSim.f2 = [l2hSim.f2; transitions(l2h).f2];
        h2lSim.f1 = [h2lSim.f1; transitions(h2l).f1];
        h2lSim.f2 = [h2lSim.f2; transitions(h2l).f2];
    end
end

tLists = {l2hRival h2lRival l2hSim h2lSim};
tLabels = {'rivalry low to high' 'rivalry high to low' 'sim low to high' 'sim high to low'};

% must have a minimum number of transitions for each run type
if size(l2hRival.f1, 1) >= minInstMod && size(h2lRival.f1, 1) >= minInstMod
    if size(l2hSim.f1, 1) >= minInstMod && size(h2lSim.f1, 1) >= minInstMod
    
        % plot averaged transitions
        for tType = 1:length(tLists)
            
            meanF1trace = mean(tLists{tType}.f1,1);
            errorF1trace = ste(tLists{tType}.f1);
            
            meanF2trace = mean(tLists{tType}.f2,1);
            errorF2trace = ste(tLists{tType}.f2);
            
            figure
            title(tLabels{tType})
            hold on
            mseb([(-transHalf+(1/512)):(1/512):transHalf],[meanF1trace; meanF2trace],[errorF1trace; errorF2trace],[],1);
            vline(0,'k','transition')
            xlabel('Time from button press (s)');
            ylabel('Amplitude');
            legend('Low Frequency', 'High Frequency');
            
            % calculate rivalry modulation index
            
            F1amp = max(meanF1trace) - min(meanF1trace);
            F2amp = max(meanF2trace) - min(meanF2trace);
            
            if tType <= 2 % rivalry runs
                rivalryF1amp(tType) = F1amp;
                rivalryF2amp(tType) = F2amp;
            else
                simF1amp(tType-2) = F1amp;
                simF2amp(tType-2) = F2amp;
            end           
        end
        
        amps = [rivalryF1amp simF1amp];
        save([indicesDir 'Amps/' parName], 'amps'); 
        
        F1index = mean(rivalryF1amp) / mean(simF1amp);
        F2index = mean(rivalryF2amp) / mean(simF2amp);
        
        modIndex = [F1index F2index];
        save([indicesDir 'modIndices/' parName], 'modIndex');
        
    else
        disp([parName 'does not have enough sim transition instances'])
        modIndex = NaN;
    end
else
    disp([parName 'does not have enough rivalry transition instances'])
    modIndex = NaN;
end



end

