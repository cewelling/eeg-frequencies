function rateFFT( parName )
%rateFFT Analyze oscillation frequency using FFT 
%   For each participant, takes smoothed RLS spectra corresponding to each FOI, subtracts low
%   frequency from high frequency, takes an FFT of the resulting time
%   course; averages over all trials

% Load parameters
analysisParams

Fs = sampRate; % hz
L = 200000; % length each trial (AFTER PADDING) in seconds

% Get EEG file info
EEGfiles = dir([eegDir parName '*']);
% ...if the participant exists...
if isempty(EEGfiles)
    return;
end
date = strtok(EEGfiles(1).date);
date = datestr(date, dateFormat);

%for iRunType = 1:length(runTypes)
for iRunType = 1 % dartRival only
    cRunType = runTypes{iRunType};
    
    ffts = [];
    
    for runIndex = runIndices
        runName = [cRunType num2str(runIndex)];
        
        %% Load RLS spectra
        
%         if exist(['rls_data/' parName '_' runName '.mat'], 'file') && isequal(analFreqs, stimFreqs)
%             load(['rls_data/' parName '_' runName '.mat'])
%         elseif exist(['rls_data/' parName '_' runName '_harmonics.mat'], 'file') && isequal(analFreqs,harFreqs)
%             load(['rls_data/' parName '_' runName '_harmonics.mat'])
%         else
%              % Get the appropriate EEG file for this run
%             EEGfile = [eegDir parName '_' runName '_' date '.bdf'];
%             
%             % skip over participants / runs that don't exist
%             if ~exist(EEGfile, 'file')
%                 continue;
%             end
%             
%             % run RLS analysis
%             [rls_data, ~] = runRLS(parName, runName, date, EEGfile);
%         end
        
        %% TEMP
         % Get the appropriate EEG file for this run
            EEGfile = [eegDir parName '_' runName '_' date '.bdf'];
            
            % skip over participants / runs that don't exist
            if ~exist(EEGfile, 'file')
                continue;
            end
        
        [rls_data, ~] = runRLSTEST(parName, runName, date, EEGfile);
        
        %%
        
        % Handle each trial in the RLS data
        for iTrial = 1:size(rls_data(1).amp, 2)
            
            traces = [];
            
            %% Subtract low frequency trace from high frequency trace (H-L)
            for iFreq = 1:2
            RLStrace = (rls_data(iFreq).amp{iTrial});
            
            % account for trials that weren't recorded properly (mistakes)
            if isempty(RLStrace)               
                continue;
            end
            
            % extra smoothing
            %HminusL = HRsmoothing(HminusL, 'gaussian', 1225, 1, 0);
            
            %% Cut off first and last second of trial
            % rivalry oscillations may be less regular on edges
            
            % trial length is different for a few trials that had a baseline period added after the start signal
            %if strcmp(parName, 'cumulus10') || strcmp(parName, 'cumulus13')
               % cutLength = sampRate*.25;
            %else
                cutLength = sampRate;
            %end
            RLStrace = RLStrace(cutLength: end - cutLength); % This does not really cut off the end since the end was padded
            
            %% Detrend by subtracting a fourth polynomial fitting
            [p,~,mu] = polyfit([1:length(RLStrace)],RLStrace, 4);
            fitting = polyval(p,[1:length(RLStrace)], [], mu);
            
            % detrended time course
            RLStrace_dt = RLStrace - fitting;
            
            %% Apply windowing function
            win = hanning(length(RLStrace_dt));
            RLStrace_win = RLStrace_dt.*win';
            
            %% Pad with zeroes to improved post-fft frequency resolution
            
            % trial length is different for a few trials that had a baseline period added after the start signal
            if strcmp(parName, 'cumulus10') || strcmp(parName, 'cumulus13')
                padLength = L/1000*512 - (28.5*sampRate - 2*cutLength);
            else
                padLength = L/1000*512 - (30*sampRate - 2*cutLength);
            end
            
            pad = zeros(1, padLength);
            RLStrace_fin = [RLStrace_win pad];
            
            
            %% Take fourier transform of the time course
            
            Y = fft(RLStrace_fin);
            
            % get one-sided fft
            P2 = angle(Y/L);
            phase = P2(1:L/2 + 1);
            phase(2:end-1) = 2*phase(2:end - 1);
            
            % store phase spectrum for this frequency
            traces = [traces; phase];
            
%             plot(freq,pow)
            
%             figure
%             plot(HminusL)
            end
            
            %% Subtract high freq phase from low freq phase
            HminusL_ph = traces(1,:) - traces(2,:);
            ffts = [ffts; HminusL_ph];
        end
        
    end
    
    %% Find the peak of the FFT plot and store the corresponding frequency
    
    % frequency axis
    freq = Fs*(0:(L/2))/L;
    
    meanfft = nanmean(ffts, 1);
    
    %     % smooth the power spectrum
%         smoothfft = HRsmoothing(meanfft, 'gaussian', 61, 1, 0);
%     
%              [~,I] = max(smoothfft);
%              oscFreq = freq(I);
    
    % Take a weighted average of the power between 0 and 0.2
%     [~, fstart] = min(abs(freq - 0));
%     [~, fend] = min(abs(freq-0.2));
%     analPiece = meanfft(fstart:fend);
    % divide each value by the sum (so that they sum to 1)
%     analPiece = analPiece/sum(analPiece);
%     freqsWted = [];
%     for i = 1:length(analPiece)
%         freqsWted = [freqsWted freq(i)*analPiece(i)];
%     end
%     oscFreq = sum(freqsWted);
    
    % Fit the peak to a quadratic equation
    %     [~, fstart] = min(abs(freq-0));
    %     [~, fend] = min(abs(freq-0.2));
    
    %     qFit = polyfit(freq(fstart:fend), meanfft(fstart:fend), 2);
    %     oscFreq = -qFit(2)/(2*qFit(1));
    
    % plot fit over data
%     figure
%     plot(freq(fstart:fend), meanfft(fstart:fend),'b.')
    %fit = qFit(1)*(freq(fstart:fend).^2) + qFit(2)*freq(fstart:fend) + qFit(3);
%     hold on
%     plot(freq(fstart:fend), fit, 'r');
    
    
    %save([indicesDir 'rateFreqs/' parName], 'oscFreq')
    
    
    %% Plot the power spectrum averaged over all trials
    
    figure
    title([parName ' fft of H-L RLS time course'])
    hold on
    %plot(freq, smoothfft);
    plot(freq, meanfft);
    xlabel('frequency (hz)')
    ylabel('phase')
    %xlim([0 2]);
%     plot(freq(fstart:fend), fit, 'r');
    
end
end

