% Runs each component of SSVEP binocular rivalry analysis of individuals
% Choose the components under "Switches" 
% Edit parameters in analysisParams.m

clearvars -except pValues kSize cutOff rValues groupNum

%add paths
addpath(genpath(['helperFunctions']))
addpath(genpath(['miscFunctions']))

%%%%% Switches: Turn on (1) or off (0) %%%%%%%%
fftAnal = 0; % plain fft, calculates SNRs
keyAnal = 0; % parses keypress data (for ampClust)
rlsAnal = 1; % rls analysis
transAnal = 0; % gets RLS traces for transitions, plots transitions per run
peakAnal = 0; % gets peaks, calculates amplitudes
oscAnal = 0; % naive analysis of oscillation frequency by FFT
ampClust = 0; % analyzes how dominant/mixed amplitudes cluster
plvAnal = 0; % analyzes the phase relationships between RLS traces
parPlot = 0; % creates participant transition and peak plots


%load parameters
paramsFlag = 'none'; % no specific parameters yet
analysisParams

% delete saved files
%delete pre-processing/preproc/*.mat
%delete rls_data/unsmoothedRLS/*.mat

% Analyze each participant
for iPar = 1:length(parNums)
    
    %%%%% Participant-specific set-up %%%%%%%%%
    
    numFormat = '%02d';
    parName = [group num2str(parNums(iPar), numFormat)];
    
    % Get EEG file info   
    EEGfiles = dir([eegDir parName '*']);
    % ...if the participant exists...
    if isempty(EEGfiles)
        continue;
    end
    date = strtok(EEGfiles(1).date);
    date = datestr(date, dateFormat);
    
    %%
    % Analyze each run (a run is a block of trials)
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
            
            %runRLS_noSmoothing(parName, runName, date, EEGfile);
            %% runFFT
            if fftAnal == 1
                paramsFlag = 'fftAnal';
                runFFT(EEGfile, runName, parName, iPar, date, paramsFlag) % add intermodulation frequencies
                %runRLS_noSmoothing(parName, runName, date, EEGfile);
            end
            % plots the FFT spectrum for each run
            % stores SNR indices (1 per stimulation frequency per run)
            
            %% Identify perceptual epochs
            if keyAnal == 1
                paramsFlag = 'none'; % no specific parameters
                [epochs, simSchedule] = getEpochs(parName, runName, date, paramsFlag);
            end
            
            %% RLS analysis
            if rlsAnal == 1
                paramsFlag = 'rlsAnal';
                [rls_data, rls_time] = runRLS(parName, runName, date, EEGfile, paramsFlag);
            end
            % Time frequency analysis (using adaptive RLS filter)
            
            %% Analyze Transitions
            if transAnal == 1
                paramsFlag = 'transAnal';
                transitions = aggTransitions(parName, iPar, runName, date, EEGfile, paramsFlag);
            end
            % Aggregate and plot transitions for each run
            
            %% Analyze Peaks
            if peakAnal == 1
                paramsFlag = 'peakAnal';
                aggPeaks(parName, runName, date, EEGfile);
            end
           
        end
    end
    
    %% Determine electrodes for analysis
    
    %% Get oscillation frequency
    if oscAnal == 1
        paramsFlag = 'oscAnal';
        rateFFT(parName, paramsFlag);
    end
    
    %% Analyze amplitudes
    if ampClust == 1
        paramsFlag = 'ampClust';
        ampClusters(parName, paramsFlag);
    end
    
    %% Create participant plots
    if parPlot == 1
        paramsFlag = 'parPlot';
        makeParPlot(parName, iPar, date, paramsFlag);
    end
    
    %% Get oscillation frequency
    if plvAnal == 1
        paramsFlag = 'oscAnal';
        computePLV(parName, cRunType, 'Shuffle', paramsFlag);
    end
    
    %% Quartile analysis
    %quartileAnalysis(parName, date);
    
    %getIMsnrs(parName);

end