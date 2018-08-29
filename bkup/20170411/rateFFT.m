function [freq, meanfft] = rateFFT( parName )
%rateFFT Analyze oscillation frequency using FFT
%   For each participant, takes smoothed RLS spectra corresponding to each FOI, subtracts low
%   frequency from high frequency, takes an FFT of the resulting time
%   course; averages over all trials

global kSize 
%global cutOff

% Load parameters
analysisParams

Fs = sampRate; % hz
L = 2000; % length of each trial (AFTER PADDING) in seconds
N = sampRate*L;


% Get EEG file info
EEGfiles = dir([eegDir parName '*']);
% ...if the participant exists...
if isempty(EEGfiles)
    return;
end
date = strtok(EEGfiles(1).date);
date = datestr(date, dateFormat);

% only take runs above min SNR
runList = {};
excluded = {};
for iRunType =  1 %1 % dartRival only
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

%for iRunType = 1:length(runTypes)
ffts = [];
f = []; % frequencies
maxes = [];

% for test-retest
oddFFTs = [];
evenFFTs = [];

if length(runList) > 1 % reject participants with only one good run
    for i = 1:length(runList)
        runName = runList{i};
        
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
            
            % is the SNR of this trial high enough?
            load(['highSNRelecs/' parName '_' runName '_trial' num2str(iTrial) '_' num2str(numElecs) 'elecs.mat']);
            if nanmean(maxSNRs(2,ismember(maxSNRs(1,:), elecs))) < minSNR
               continue;
            end
            
            %% Subtract low frequency trace from high frequency trace (H-L)
            
            RLSlow = (rls_data(1).amp{iTrial});
            RLShigh = (rls_data(2).amp{iTrial});
            
            % account for trials that weren't recorded properly (mistakes)
            if isempty(RLSlow)
                continue;
            end
            
            HminusL = RLShigh - RLSlow;
            
            % Create sine wave for testing
                     %frequency = 0.50; %hz
                     %HminusL = 2*sin(frequency*2*pi*(1:1/512:30));
                     %HminusL = [HminusL(1:4000) zeros(1,6000) HminusL(10000:14000)];
            
            % extra smoothing
%             if kSize ~= 1
%                 HminusL = HRsmoothing(HminusL, 'gaussian', kSize, 1, 0);
%             end
            %HminusL = HRsmoothing(HminusL, 'gaussian', 201, 1, 0);
            %HminusL = HRsmoothing(HminusL, 'gaussian', smooth_ksize, smooth_sd, 0);
            
            %% Filtering
             %if cutOff ~= 0
%                 Hd = designfilt('lowpassfir','FilterOrder',5000,'CutoffFrequency',0.25, ...
%                     'DesignMethod','window','Window','hamming','SampleRate',Fs);
%                 % compensate for filter delay
%                 D = mean(grpdelay(Hd));
%                 HminusL = filter(Hd, [HminusL zeros(1,D)]);
%                 HminusL = HminusL(D+1:end);
             %end
            
            %% Cut off first and last second of trial
            % rivalry oscillations may be less regular on edges
            
            % trial length is different for a few trials that had a baseline period added after the start signal
            %if strcmp(parName, 'cumulus10') || strcmp(parName, 'cumulus13')
            % cutLength = sampRate*.25;
            %else
            cutLength = sampRate;
            %end
            HminusL = HminusL(cutLength: end - cutLength); % This does not really cut off the end since the end was padded
            
            %% Detrend by subtracting a fourth polynomial fitting
            [p,~,mu] = polyfit([1:length(HminusL)],HminusL, 4);
            fitting = polyval(p,[1:length(HminusL)], [], mu);
            %fitting = mean(HminusL); % demean
            
            % detrended time course
            HminusL_dt = HminusL - fitting;
            
            %% Apply windowing function
            %win = hanning(length(HminusL_dt));
            %HminusL_win = HminusL_dt.*win';
            
            
            %% Take fourier transform of the time course
            
            %         %Y = fft(HminusL,N);
            %         Y = fft(HminusL_dt,N);
            %         %figure
            %         %plot(HminusL)
            %         % get one-sided fft
            %         P2 = abs(Y/N);
            %         pow = P2(1:N/2 + 1);
            %         pow(2:end-1) = 2*pow(2:end - 1);
            %
            %         % store fft spectra for this type of transition
            %         ffts = [ffts; pow];
            %
            %         %plot(freq,pow)
            %
            %         %figure
            %         %plot(HminusL_fin)
            
            %% Use chronus tool box to get multi-taper spectrum
            
            % Ensure that path is in correct order
            %path(path, genpath('/Users/acspiegel/Documents/MATLAB/Toolboxes/FieldTrip'));

            tw = 3; %2;
            %         if length(HminusL_dt) < sampRate * (trialDur - 2) % runs with baseline period after start signal
            %             tw = (trialDur - 1.5 - 2) * bwidth;
            %         else
            %             tw = (trialDur - 2) * bwidth;
            %         end
            params.tapers = [tw tw*2 - 1];
            %params.tapers = [30 59];
            params.Fs = 512;
            params.pad = 4;
            %params.fpass = [0 0.5];
            
            [S, f] = mtspectrumc(HminusL_dt, params);
            %         figure
            %         plot(f,S)
            %         xlim([0 3])
            %         title(['tw = ' num2str(tw)])
            %ffts = [ffts; S'/sum(S)];
            ffts = [ffts; S'];
            if mod(iTrial,2) > 0
                oddFFTs = [oddFFTs; S'];
            else
                evenFFTs = [evenFFTs; S'];
            end
            [~,I] = max(S);
            maxes = [maxes; f(I)];
            
        end
    end
end



%% Find the peak of the FFT plot and store the corresponding frequency

if ~isempty(ffts)
    
    % frequency axis
    %freq = Fs*(0:(N/2))/N; %freq = 0:1/(L/1000):Fs/2;
    freq = f;
    
    meanfft = nanmean(ffts, 1);
    oMeanFFT = nanmean(oddFFTs, 1);
    eMeanFFT = nanmean(evenFFTs, 1);
    oNum = size(oddFFTs, 1);
    eNum = size(evenFFTs, 1);
    
%     % calculate the cumulative distribution function
%     CDF = cumsum(meanfft);
%     CDF = CDF/max(CDF);
%     % find the index of the half max
%     [~, hMaxI] = min(abs(.5 - CDF));
%     oscFreq = freq(hMaxI);
%     oscAmp = meanfft(hMaxI);
    
    %     % smooth the power spectrum
    %smoothfft = HRsmoothing(meanfft, 'gaussian', 401, 1, 0);
    
    %[~,I] = max(smoothfft);
    [A,I] = max(meanfft);
    oscFreq = freq(I);
    %oscFreq = mean(maxes);
    oscAmp = A;
    
    [~, oddI] = max(oMeanFFT);
    oscFreq_odd = freq(oddI);
    [~, evenI] = max(eMeanFFT);
    oscFreq_even = freq(evenI);
    
    %% Plot the power spectrum averaged over all trials
    
%         figure
%         subplot(2,1,1)
%         title([parName ' fft of H-L RLS time course'])
%         hold on
%         plot(freq, meanfft);
%         %plot(freq, ffts');
%         %plot(freq, smoothfft, 'LineWidth',3);
%         %legend('unsmoothed','smoothed')
%         vline(oscFreq)
%         xlabel('frequency (hz)')
%         ylabel('power')
%         xlim([0 2]);
%         %plot(freq(fstart:fend), fit, 'r');
%         %ylim([0 10])
%         
% %         subplot(3,1,2)
% %         plot(freq,CDF);
% %         xlabel('frequency (hz)')
% %         xlim([0 2])
% %         vline(oscFreq)
%         
%         subplot(2,1,2)
%         plot(HminusL_dt)
%         %ylim([-3 3])
    
else
    oscFreq = NaN;
    oscAmp = NaN;
    freq = [];
    meanfft = [];
    
    oscFreq_odd = NaN;
    oscFreq_even = NaN;
    oNum = NaN;
    eNum = NaN;
end

% Take a weighted average of the power between 0 and 0.2
%         [~, fstart] = min(abs(freq - 0));
%         [~, fend] = min(abs(freq-0.2));
%         analPiece = meanfft(fstart:fend);
% divide each value by the sum (so that they sum to 1)
%         analPiece = analPiece/sum(analPiece);
%         freqsWted = [];
%         for i = 1:length(analPiece)
%             freqsWted = [freqsWted freq(i)*analPiece(i)];
%         end
%         oscFreq = sum(freqsWted);

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

if iRunType == 2
    save([indicesDir 'rateFreqs/' parName '_SIM'], 'oscFreq')
    save([indicesDir 'oscAmps/' parName '_SIM'], 'oscAmp')
else
    save([indicesDir 'rateFreqs/' parName], 'oscFreq')
    save([indicesDir 'oscAmps/' parName], 'oscAmp')
    
    % for test-retest
    save([indicesDir 'rateFreqs/' parName '_split'], 'oscFreq_odd', 'oscFreq_even', 'oNum', 'eNum');
end

oscFreq

end

