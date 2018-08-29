function [freq, meanfft] = rateFFT( parName )
%rateFFT Analyze oscillation frequency using FFT
%   For each participant, takes smoothed RLS spectra corresponding to each FOI, subtracts low
%   frequency from high frequency, takes an FFT of the resulting time
%   course; averages over all trials

%global kSize
global cutOff

clearvars -except kSize parName cutOff

% delete all existing oscillation frequencies in folder for this
% participant
delete(['indices/rateFreqs/' parName '*.mat'])

% Load parameters
analysisParams

Fs = sampRate; % hz
L = 2000; % length of each trial (AFTER PADDING) in seconds
N = sampRate*L;


% Get EEG file info
EEGfiles = dir([eegDir parName '*']);
% ...if the participant exists...
if isempty(EEGfiles)
    freq = [];
    meanfft = [];
    return;
end
date = strtok(EEGfiles(1).date);
date = datestr(date, dateFormat);

% only take trials above min SNR
excluded = zeros(3,numTrials);

for iRunType =  1; %1 % dartRival only
    cRunType = runTypes{iRunType};
    
    ffts = [];
    f = []; % frequencies
    maxes = [];
    
    % for test-retest
    oddFFTs = [];
    evenFFTs = [];
    splitGroup = 1;
    
    trialGroup = 1; % initialize trial group, for test-retest reliability purposes
    
    for runIndex = runIndices
        runName = [cRunType num2str(runIndex)];
        
        %% Load RLS spectra
        
        % Get the appropriate EEG file for this run
        EEGfile = [eegDir parName '_' runName '_' date '.bdf'];
        
        % skip over participants / runs that don't exist
        if ~exist(EEGfile, 'file')
            continue;
        end
        
        [rls_data, rls_time] = runRLS_noSmoothing(parName, runName, date, EEGfile);
        
        %%
        
        % Handle each trial in the RLS data
        for iTrial = 1:size(rls_data(1).amp, 2)
            
            %% Subtract low frequency trace from high frequency trace (H-L)
            
            RLSlow = (rls_data(1).amp{iTrial});
            RLShigh = (rls_data(2).amp{iTrial});
            
            % account for trials that weren't recorded properly (mistakes)
            if isempty(RLSlow)
                continue;
            end
            
            HminusL = RLShigh - RLSlow;
            
            %% Is the SNR of this trial high enough?
            load(['pre-processing/highSNRelecs/' parName '_' runName '_trial' num2str(iTrial) '_' electrodeSet.name '.mat']);
            %if nanmean(maxSNRs(2,ismember(maxSNRs(1,:), elecs))) < minSNR
%             if nanmean(lFreqSNRs(2,ismember(maxSNRs(1,:), elecs))) < 2 || ...
%                     nanmean(hFreqSNRs(2,ismember(maxSNRs(1,:), elecs))) < 2
%            if nanmean(maxSNRs(2,maxSNRs(1,:)==snrElecs)) < minSNR
             if nanmean(lFreqSNRs(2,lFreqSNRs(1,:)==snrElecs)) < minSNR || ...
                     nanmean(hFreqSNRs(2,hFreqSNRs(1,:)==snrElecs)) < minSNR

                
                % record excluded trial
                excluded(runIndex,iTrial) = 1;
                continue;
            end
            
            %% Choose which trace to use
            
            %             if nanmean(lFreqSNRs(2,ismember(maxSNRs(1,:), elecs))) > hFreqSNRs(2,ismember(maxSNRs(1,:), elecs))
            %                 HminusL = RLSlow;
            %             else
            %                 HminusL = RLShigh;
            %             end
            
            % Create sine wave for testing
            %                      frequency = 0.8; %hz
            %                      HminusL = 2*sin(frequency*2*pi*(1:1/512:30));
            %                      HminusL = [HminusL(1:4000) zeros(1,6000) HminusL(10000:14000)];
            
            % extra smoothing
            %             if kSize ~= 1
            %                 HminusL = HRsmoothing(HminusL, 'gaussian', kSize, 1, 0);
            %             end
            %HminusL = HRsmoothing(HminusL, 'gaussian', 201, 1, 0);
            %HminusL = HRsmoothing(HminusL, 'gaussian', smooth_ksize, smooth_sd, 0);
            
            %% Filtering
            %                         if cutOff ~= 0
            %                             Hd = designfilt('lowpassfir','FilterOrder',5000,'CutoffFrequency',0.7, ...
            %                                 'DesignMethod','window','Window','hamming','SampleRate',Fs);
            %                             % compensate for filter delay
            %                             D = mean(grpdelay(Hd));
            %                             HminusL = filter(Hd, [HminusL zeros(1,D)]);
            %                             HminusL = HminusL(D+1:end);
            %                         end
            
            %% Cut off first and last second of trial
            % rivalry oscillations may be less regular on edges
            
            % trial length is different for a few trials that had a baseline period added after the start signal
            %if strcmp(parName, 'cumulus10') || strcmp(parName, 'cumulus13')
            % cutLength = sampRate*.25;
            %else
            cutLength = sampRate;
            %end
            HminusL = HminusL(cutLength: end - sampRate); % This does not really cut off the end since the end was padded
            
            %% Detrend by subtracting a fourth polynomial fitting
            %[p,~,mu] = polyfit([1:length(HminusL)],HminusL, 4);
            %fitting = polyval(p,[1:length(HminusL)], [], mu);
            fitting = mean(HminusL); % demean
            
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
            
            tw = 2; %3: 20170525; %2;
            %         if length(HminusL_dt) < sampRate * (trialDur - 2) % runs with baseline period after start signal
            %             tw = (trialDur - 1.5 - 2) * bwidth;
            %         else
            %             tw = (trialDur - 2) * bwidth;
            %         end
            params.tapers = [tw tw*2 - 1];
            %params.tapers = [30 59];
            params.Fs = sampRate;
            params.pad = 4;
            %params.fpass = [0 0.5];
            
            try
                [S, f] = mtspectrumc(HminusL_dt, params);
            catch % fieldtrip cannot be at the top of the path for this function to work
                path(path, genpath('/Users/acspiegel/Documents/MATLAB/Toolboxes/FieldTrip'));
                [S, f] = mtspectrumc(HminusL_dt, params);
            end
            
            %         figure
            %         plot(f,S)
            %         xlim([0 3])
            %         title(['tw = ' num2str(tw)])
            %ffts = [ffts; S'/sum(S)];
            ffts = [ffts; S'];
            
%             % cumulus 29 doesn't have any good odd trials
%             if strcmp(parName, 'cumulus29')
%                 if ((iTrial == 2 || iTrial == 6) && runIndex == 1) || ...
%                         (iTrial == 4 && runIndex == 3)
%                     oddFFTs = [oddFFTs; S'];
%                 else
%                     evenFFTs = [evenFFTs; S'];
%                 end
%             else
%                 if mod(iTrial,2) > 0
%                     oddFFTs = [oddFFTs; S'];
%                 else
%                     evenFFTs = [evenFFTs; S'];
%                 end
%             end

            if splitGroup == 1
                 oddFFTs = [oddFFTs; S'];
                 splitGroup = 2;
            else
                 evenFFTs = [evenFFTs; S'];
                 splitGroup = 1;
            end
            
            [~,I] = max(S);
            maxes = [maxes; f(I)];
            
        end
    end
    
    % save excluded trials (for rivalry trials only)
    if iRunType == 1
        save(['pre-processing/badTrials/' parName], 'excluded')
    end
    
    %% Find the peak of the FFT plot and store the corresponding frequency
    
    if size(ffts) >= 6
        
        % frequency axis
        %freq = Fs*(0:(N/2))/N; %freq = 0:1/(L/1000):Fs/2;
        freq = f;
        
        meanfft = nanmean(ffts, 1);
        oMeanFFT = nanmean(oddFFTs, 1);
        eMeanFFT = nanmean(evenFFTs, 1);
        oNum = size(oddFFTs, 1);
        eNum = size(evenFFTs, 1);
        
        %     % calculate the cumulative distribution function
        CDF = cumsum(meanfft);
        CDF = CDF/max(CDF);
        % find the index of the half max
        [~, hMaxI] = min(abs(.5 - CDF));
        oscFreq = freq(hMaxI);
        oscAmp = meanfft(hMaxI);
        
        %     % smooth the power spectrum
        smoothfft = HRsmoothing(meanfft, 'gaussian', 31, 1, 0);
        
        %%%[A,I] = max(meanfft);
        %[A,I] = max(smoothfft);
        %%%oscFreq = freq(I);
        %oscFreq = mean(maxes);
        %%%oscAmp = A;
        
        %% For test-retest reliability
        
        % calculate the cumulative distribution function
        oCDF = cumsum(oMeanFFT);
        oCDF = oCDF/max(oCDF);
        % find the index of the half max
        [~, ohMaxI] = min(abs(.5 - oCDF));
        oscFreq_odd = freq(ohMaxI);
        
        eCDF = cumsum(eMeanFFT);
        eCDF = eCDF/max(eCDF);
        % find the index of the half max
        [~, ehMaxI] = min(abs(.5 - eCDF));
        oscFreq_even = freq(ehMaxI);
        
        %% Save the power spectrum and CDF
        if iRunType == 2
            save(['oscFFTs/' parName '_SIM'], 'freq', 'meanfft', 'CDF')
        else
            save(['oscFFTs/' parName], 'freq', 'meanfft', 'CDF')
        end
        
        %% Plot the power spectrum averaged over all trials
        
        figure
        %subplot(2,1,2)
        if iRunType == 2
            title([parName ' fft of H-L SIM RLS time course    .'])
        else
            title([parName ' fft of H-L RLS time course'])
        end
        hold on
        plot(freq, meanfft);
        %plot(freq, ffts');
        %plot(freq, smoothfft, 'LineWidth',3);
        plot(freq, CDF);
        %legend('unsmoothed','smoothed')
        %vline(oscFreq)
        %plot([oscFreq oscFreq],[0 meanfft(I)],'r-')
        plot([oscFreq oscFreq],[0 meanfft(hMaxI)],'r-')
        xlabel('Frequency (Hz)')
        ylabel('Power')
        xlim([0 2]);
        %plot(freq(fstart:fend), fit, 'r');
        %ylim([0 10])
        
        %         subplot(3,1,2)
        %         plot(freq,CDF);
        %         xlabel('frequency (hz)')
        %         xlim([0 2])
        %         vline(oscFreq)
        
        %subplot(2,1,1)
        figure
        plot(rls_time(cutLength:end-cutLength),HminusL_dt)
        box off
        xlabel('time (s)')
        ylabel('amplitude')
        legend('8.5 Hz minus 5.67 Hz', 'Location', 'NorthEast');
        legend boxoff
        ylim([-3 3])
        
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
end


