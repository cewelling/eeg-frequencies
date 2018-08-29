% Runs each component of SSVEP binocular rivalry analysis of invdividuals
% Choose those components under "Switches" 
% Edit parameters in analysisParams.m

clearvars

%%%%% Switches: Turn on (1) or off (0) %%%%%%%%
fftAnal = 1; %plain fft, calculates SNRs
keyAnal = 0; %parses keypress data
durAnal = 0; %calculates epoch durations, suppression indices
rlsAnal = 0; %rls analysis
transAnal = 0; %plots transitions
peakAnal = 0; %gets peaks, calculates amplitudes
transCount = 0; %counts transitions
modAnal = 0; %calculate rivalry modulation indices
getElecs = 0; %determine analysis electrodes based on SNR
parPlot = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES FOR modAnal params: (edit these in analysisParams.m)
% 1. Ensure that you have selected both rivalry and sim runs to compare
% 2. Ensure that you have selected the same electrode set for SNR and RLS analysis
% 3. Choose just one noise window half for SNR analysis
% 4. RUN STRATUS134 separately (only uses runNums 1-5 for dartboards)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diffAnal = 0; %difference scores between two EEG freqs
alignAnal = 0; %aligns keypress with EEG


%load parameters
analysisParams

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
            
            %% runFFT
            if fftAnal == 1
                runFFT(EEGfile, runName, parName, iPar, date) % add intermodulation frequencies
            end
            % plots the FFT spectrum for each run
            % stores SNR indices (1 per stimulation frequency per run)
            
            %% Identify perceptual epochs
            if keyAnal == 1
                [epochs, simSchedule] = getEpochs(parName, runName, date);
            end
            
            %% Get epoch durations
            if durAnal == 1
                getEpochDurations(parName, runName, date)
            end
            % store epoch duration indices (1 per epoch type (i.e. high freq dom, low freq dom, mixed) per run)
            % store suppression indices--fraction of time spent in a dominant state
            
            %% RLS analysis
            if rlsAnal == 1
                [rls_data, rls_time] = runRLS(parName, runName, date, EEGfile);
            end
            % Time frequency analysis (using adaptive RLS filter)
            
            %% Analyze Transitions
            if transAnal == 1
                transitions = aggTransitions(parName, runName, date, EEGfile);
            end
            % Aggregate and plot transitions for each run
            
            %% Analyze Peaks
            if peakAnal == 1
                aggPeaks(parName, runName, date, EEGfile);
            end
            
            %% Calculate transition rate
            if transCount == 1
                tCount = getTransitions(parName, runName, date);
            end
        end
    end
    
    %% Compare modulation
    if modAnal == 1
        getModIndex(parName, date, EEGfiles, iPar)
    end
    % Calculate a rivalry modulation index comparing response to rivalry to
    % response to stimulus alternation (sim trials) for each participant
    
    %% Determine electrodes for analysis
    % electrodes used are the same for all of each participant's dart runs,
    % all of each participant's marz runs
    if getElecs == 1
        setElecs(parName);
    end
    
    %% Create participant plots
    if parPlot == 1
        makeParPlot(parName, date, EEGfiles, iPar);
    end

end

% Group analysis