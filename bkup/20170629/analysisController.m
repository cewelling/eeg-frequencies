% Runs each component of SSVEP binocular rivalry analysis of individuals
% Choose the components under "Switches" 
% Edit parameters in analysisParams.m

clearvars -except pValues kSize cutOff rValues groupNum

%%%%% Switches: Turn on (1) or off (0) %%%%%%%%
fftAnal = 0; % plain fft, calculates SNRs
getElecs = 0; % determines analysis electrodes based on SNR
keyAnal = 0; % parses keypress data (for ampClust)
rlsAnal = 0; % rls analysis
transAnal = 0; % plots transitions
peakAnal = 0; % gets peaks, calculates amplitudes
oscAnal = 0; % naive analysis of oscillation frequency by FFT
ampClust = 0; % analyzes how dominant/mixed amplitudes cluster
parPlot = 1; % creates participant transition and peak plots


%load parameters
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
                runFFT(EEGfile, runName, parName, iPar, date) % add intermodulation frequencies
                %runRLS_noSmoothing(parName, runName, date, EEGfile);
            end
            % plots the FFT spectrum for each run
            % stores SNR indices (1 per stimulation frequency per run)
            
            %% Identify perceptual epochs
            if keyAnal == 1
                [epochs, simSchedule] = getEpochs(parName, runName, date);
            end
            
            %% RLS analysis
            if rlsAnal == 1
                [rls_data, rls_time] = runRLS(parName, runName, date, EEGfile);
            end
            % Time frequency analysis (using adaptive RLS filter)
            
            %% Analyze Transitions
            if transAnal == 1
                transitions = aggTransitions(parName, iPar, runName, date, EEGfile);
            end
            % Aggregate and plot transitions for each run
            
            %% Analyze Peaks
            if peakAnal == 1
                aggPeaks(parName, runName, date, EEGfile);
            end
        end
    end
    
    %% Determine electrodes for analysis
    % electrodes used are the same for all of each participant's dart runs, all of each participant's marz runs
    if getElecs == 1
        setElecs(parName);
    end
    
    %% Get oscillation frequency
    if oscAnal == 1
        rateFFT(parName);
    end
    
    %% Analyze amplitudes
    if ampClust == 1
        ampClusters(parName);
    end
    
    %% Create participant plots
    if parPlot == 1
        makeParPlot(parName, iPar, date);
    end
    
    %% Quartile analysis
    %quartileAnalysis(parName, date);
    
    %getIMsnrs(parName);

end