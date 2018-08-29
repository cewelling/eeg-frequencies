function [ peaks ] = aggPeaks( parName, runName, date, EEGfile)
% aggPeaks(parName, runName, date, EEGfile) finds and saves peaks in the 
% pLists folder. Also saves chunks of RLS data corresponding to those peaks
% in the pListsRLS folder.
%
% [peaks] = aggPeaks(parName, runName, date, EEGfile) returns a
% struct with entries corresponding low and high frequency peaks, 
% respectively. 
%
% peaks struct fields: 
% type: lowFreq or highFreq
% f1: low frequency traces
% f2: high frequency traces
% durs: epoch durations
% snrs: snrs of the trials in which the peaks occurred
%
% Plots RLS traces with identified peaks overlaid. Also plots average of 
% peaks for the run designated. Note that peaks are not excluded on the 
% basis of SNR at this point.
%
% Called from: analysisController.m, makeParPlot.m
% Dependencies: analysisParams.m, runRLS.m, getEpochs.m, findPeak.m,
% (peakFinder.m)

%% Set-up

% load parameters
analysisParams

%% Prepare to label transitions with trial-by-trial SNR values (for exclusion later)

% Get stored analysis electrodes 
 if strfind(runName, 'dart')
     load(['pre-processing/highSNRelecs/' parName '_darts_' num2str(commonElecs) 'of' num2str(numElecs) 'elecs']);
 else
     load(['pre-processing/highSNRelecs/' parName '_marz_' num2str(commonElecs) 'of' num2str(numElecs) 'elecs']);
 end
 
% load SNRs for each trial in this run
snrVect = zeros(1, numTrials);
for iTrial = 1:numTrials
    if exist(['pre-processing/highSNRelecs/' parName '_' runName '_trial' num2str(iTrial) '_' num2str(numElecs) 'elecs.mat'],'file')
        load(['pre-processing/highSNRelecs/' parName '_' runName '_trial' num2str(iTrial) '_' num2str(numElecs) 'elecs.mat']);
        snrVect(iTrial) = maxSNRs(2, ismember(maxSNRs(1,:),elecs));
    else
        snrVect(iTrial) = NaN;
    end
end

%% load smoothed RLS time amplitude data
if ~exist(['rls_data/smoothedRLS/' parName '_' runName '.mat'], 'file')
    [rls_data, rls_time] = runRLS(parName, runName, date, EEGfile);
else
    load(['rls_data/smoothedRLS/' parName '_' runName '.mat'])
end

%% get list of participant's perceptual epochs for this run
    
    [epochs, ~] = getEpochs(parName, runName, date);
    % Epoch structure reminder: 1: Trial; 2: Start; 3: End; 4: button pressed; 5: dominant frequency;

%% Find and store dominant state peaks 

% Only find peaks if they stored peak lists don't already exist
if altPeak && ~exist(['peaks/altpLists/' parName '_' runName '.mat'], 'file') || ...
        ~altPeak && ~exist(['peaks/pLists/' parName '_' runName '.mat'], 'file')
    
    % Preallocate space    
    fastIndices = [];
    slowIndices = [];
    fastTrials = [];
    slowTrials = [];
    fastDurs = [];
    slowDurs = [];
    
    if blindPeak == 1
        fastTpoints = [];
        slowTpoints = [];
        fastTPtrials = [];
        slowTPtrials = [];
    end
    
    % Find peaks based on perceptual report 
   if blindPeak == 0   
       
        % Find a reported dominant percept
        epochStart = NaN;
        for iEpoch = 2:(size(epochs,1) - 1) % in first and last epoch, peaks might be cut off
            if epochs(iEpoch, 5) ~= 0
                
                % gap after epoch is long enough to signify end of epoch, or this is the last epoch in this trial
                % (combines epochs if participant quickly lifted a finger)
                if epochs(iEpoch + 1, 2) - epochs(iEpoch, 3) > gapMax || epochs(iEpoch + 1, 1) ~= epochs(iEpoch, 1)...
                        || epochs(iEpoch + 1, 5) ~= epochs(iEpoch, 5)
                    
                    if isnan(epochStart)
                        peakEpoch = epochs(iEpoch, :);
                    else
                        peakEpoch = [epochs(iEpoch, 1) epochStart epochs(iEpoch, 3:5)];
                    end
                    
                    % Step from the edges of the dom percept & find its peak
                    peakIndex = findPeak(peakEpoch, rls_data, rls_time);
                    
                    % Store peaks, trials in which they appear, and their durations
                    if ~isnan(peakIndex)
                        if epochs(iEpoch, 5) == 1 % fast freq dominant
                            fastIndices = [fastIndices peakIndex];
                            fastTrials = [fastTrials peakEpoch(1)];
                            fastDurs = [fastDurs peakEpoch(3) - peakEpoch(2)];
                        else % slow freq dominant
                            slowIndices = [slowIndices peakIndex];
                            slowTrials = [slowTrials peakEpoch(1)];
                            slowDurs = [slowDurs peakEpoch(3) - peakEpoch(2)];
                        end
                    end
                    
                    epochStart = NaN; % reset epoch start
                    
                else
                    % set start of epoch (unless start has already been set)
                    if isnan(epochStart)
                        epochStart = epochs(iEpoch, 2);
                    end
                end
            end
        end
        
    % alternate "blind" peak-picking (haoran's scripts, CURRENTLY UNFINISHED)
    else
        for iTrial = 1:length(rls_data(1).amp)
            % find transition points
            lowFderiv = diff(rls_data(1).amp{iTrial}) * sampRate;
            highFderiv = diff(rls_data(2).amp{iTrial}) * sampRate;
            trialSlowTP = peakfinder(lowFderiv, [], [], [], false);
            trialFastTP = peakfinder(highFderiv, [], [], [], false);
            slowTpoints = [slowTpoints trialSlowTP];
            fastTpoints = [fastTpoints trialFastTP];
            slowTPtrials = [slowTPtrials iTrial*ones(1,length(trialSlowTP))];
            fastTPtrials = [fastTPtrials iTrial*ones(1,length(trialFastTP))];
            
            % find peaks
            trialSlowPeaks = peakfinder(rls_data(1).amp{iTrial}, [], [], [], false);
            trialFastPeaks = peakfinder(rls_data(2).amp{iTrial}, [], [], [], false);
            
            % limit peak width (CONTINUE HERE)
            trialSlowPeaks = [];
            
            slowIndices = [slowIndices trialSlowPeaks];
            fastIndices = [fastIndices trialFastPeaks];
            slowTrials = [slowTrials iTrial*ones(1,length(trialSlowPeaks))];
            fastTrials = [fastTrials iTrial*ones(1,length(trialFastPeaks))];
        end
    end
    peakIndices = {slowIndices fastIndices};
    peakTrials = {slowTrials fastTrials};
    peakLengths = {slowDurs fastDurs};
    peakSNRs = {snrVect(slowTrials(:)) snrVect(fastTrials(:))};
    
    % save peak identification (based on participant report)
    if altPeak
        save(['peaks/altpLists/' parName '_' runName], 'peakIndices', 'peakTrials', 'peakLengths')
    else
        save(['peaks/pLists/' parName '_' runName], 'peakIndices', 'peakTrials', 'peakLengths')
    end
    
elseif altPeak && exist(['peaks/altpLists/' parName '_' runName '.mat'], 'file')
    load(['peaks/altpLists/' parName '_' runName '.mat'])
else
    load(['peaks/pLists/' parName '_' runName '.mat'])
end

%% Plot each trial with peaks marked

if strcmp(PkPlotOrNot, 'yes')
    for iTrial = 1:size(rls_data(1).amp, 2)
        
        % plot RLS amplitudes
        lowBand = rls_data(1).amp{iTrial};
        highBand = rls_data(2).amp{iTrial};
        
        figure
        title(['Peaks: ', parName ',' runName ', trial ' num2str(iTrial)])
        hold on
        plot(rls_time,lowBand,'b','linewidth',2)
        plot(rls_time,highBand,'r','linewidth',2)
        
        % find highest value that will be plotted
        checkMax = [max(lowBand) max(highBand)];
        shadeLim = ceil(max(checkMax)); % round up
        
        % overlay epochs
        trialEpochs = epochs(epochs(:,1) == iTrial,:);
        
        for iEpoch = 1:size(trialEpochs,1)
            
            horiz = [trialEpochs(iEpoch,2) trialEpochs(iEpoch,3)];
            vert = [shadeLim shadeLim];
            
            if trialEpochs(iEpoch,5) > 0.5
                area(horiz,vert,'facecolor','r','facealpha',0.15,'edgealpha',0.1);
            elseif trialEpochs(iEpoch,5) < -0.5
                area(horiz,vert,'facecolor','b','facealpha',0.15,'edgealpha',0.1);
            else
                area(horiz,vert,'facecolor','y','facealpha',0.15,'edgealpha',0.1);
            end
        end
        
        % overlay peaks
        trialFPiis = find(fastTrials == iTrial);
        trialFPis = fastIndices(trialFPiis);
        plot(rls_time(trialFPis), highBand(trialFPis), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
        trialSPiis = find(slowTrials == iTrial);
        trialSPis = slowIndices(trialSPiis);
        plot(rls_time(trialSPis), lowBand(trialSPis), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    end
end

%% Set the baseline of the RLS data

nRLS_data = [];

% If selected, get baseline using average across adjacent frequencies 
if strcmp(normType, 'freqWin') && exist([blineDir parName '_' runName '_' blineLoc '.mat'], 'file')
    load([blineDir parName '_' runName '_' blineLoc '.mat']);
elseif strcmp(normType, 'freqWin')
    baselines = getBaseline(parName, runName);
end

for iTrial = 1:size(rls_data(1).amp, 2)
    for iFreq = 1:length(rls_data);
        
        RLSamp = (rls_data(iFreq).amp{iTrial})';
        
        % account for trials that weren't recorded properly (mistakes)
        if isempty(RLSamp)
            nRLSamp = [];
            continue;
        end
        
        % take normalizing/standardizing step
        if strcmp(normType, 'mean')
            nRLSamp = RLSamp - mean(RLSamp);
        elseif strcmp(normType, 'norm')
            nRLSamp = 2*(RLSamp - min(RLSamp))/(max(RLSamp) - min(RLSamp)) - 1;
        elseif strcmp(normType, 'z')
            nRLSamp = zscore(RLSamp);
        elseif strcmp(normType, 'none')
            nRLSamp = RLSamp;
        elseif strcmp(normType, 'noStimBase')
            nRLSamp = RLSamp - baselines(iFreq);
        elseif strcmp(normType, 'freqWin')
            nRLSamp = RLSamp ./ nanmean(baselines{iTrial}{iFreq})';
        end
        nRLS_data(iFreq).amp{iTrial} = nRLSamp';
    end
    %     figure
    %     hold on
    %     lowBand = nRLS_data(1).amp{iTrial};
    %         highBand = nRLS_data(2).amp{iTrial};
    %     plot(rls_time,lowBand,'b','linewidth',2)
    %         plot(rls_time,highBand,'r','linewidth',2)
end

%% Average peaks

halfLength = sampRate*peakHalf;
fullLength = halfLength*2+1; % twice the half length plus the transition point
peakTypes = {'lowFreq', 'highFreq'};

for iPeakFreq = 1:length(peakIndices) % for peaks in both stim freqs
    freqPeaks = peakIndices{iPeakFreq};
    trials = peakTrials{iPeakFreq}; %list of trials containing this transition
    
    f1Epochs = []; % low stim freq
    f2Epochs = []; % high stim freq
    
    for iPeak = 1:length(trials)
        
        peakIndex = freqPeaks(iPeak);
        
        % get RLS data from appropriate trial
        RLSf1 = nRLS_data(1).amp{trials(iPeak)};
        RLSf2 = nRLS_data(2).amp{trials(iPeak)};
        
        % extract RLS amplitudes around peak
        
        % epoch is cut off by beginning, do not include it (to avoid discontinuities in transition traces) 
        if (peakIndex - halfLength) < 1
            epochRLSf1 = nan(1, 2*halfLength + 1);
            epochRLSf2 = nan(1, 2*halfLength + 1);
            
%             % pad beginning with NaNs
%             epochRLSf1 = RLSf1(1:peakIndex + halfLength);
%             epochRLSf2 = RLSf2(1:peakIndex + halfLength);
%             padSize = fullLength - length(epochRLSf1);
%             
%             epochRLSf1 = [nan(1,padSize), epochRLSf1];
%             epochRLSf2 = [nan(1,padSize), epochRLSf2];
            
        % epoch is cut off by beginning, do not include it (to avoid discontinuities in transition traces) 
        elseif (peakIndex + halfLength) > length(RLSf1) - (cutoffWindow / 1000 * 512); % remember end is padded
            epochRLSf1 = nan(1, 2*halfLength + 1);
            epochRLSf2 = nan(1, 2*halfLength + 1);
            
%             % pad beginning with NaNs
%             epochRLSf1 = RLSf1(peakIndex - halfLength : length(RLSf1) - (cutoffWindow / 1000 * 512));
%             epochRLSf2 = RLSf2(peakIndex - halfLength : length(RLSf2) - (cutoffWindow / 1000 * 512));
%             padSize = fullLength - length(epochRLSf1);
%             
%             epochRLSf1 = [epochRLSf1,nan(1,padSize)];
%             epochRLSf2 = [epochRLSf2,nan(1,padSize)];
            
        % epoch is not cut off
        else
            epochRLSf1 = RLSf1(peakIndex-halfLength:peakIndex + halfLength);
            epochRLSf2 = RLSf2(peakIndex-halfLength:peakIndex + halfLength);
        end
        
        f1Epochs = [f1Epochs; epochRLSf1];
        f2Epochs = [f2Epochs; epochRLSf2];
        
    end
    
    % Store peaks with properties 
    peaks(iPeakFreq).type = peakTypes{iPeakFreq};
    peaks(iPeakFreq).f1 = f1Epochs;
    peaks(iPeakFreq).f2 = f2Epochs;
    peaks(iPeakFreq).durs = peakLengths{iPeakFreq}';
    peaks(iPeakFreq).snrs = peakSNRs{iPeakFreq}';
    
end

% save peaks
if altPeak
    save(['peaks/altpListsRLS/' parName '_' runName '_' normType], 'peaks');
else
    if analFreqs == stimFreqs
        save(['peaks/pListsRLS/' parName '_' runName '_' normType], 'peaks');
    elseif analFreqs == harFreqs
        save(['peaks/pListsRLS/' parName '_' runName '_' normType '_harmonics'], 'peaks');
    end
end

%% Plot averaged peaks

if strcmp(avgPksPlotOrNot, 'yes')
    for iPeakFreq = 1:size(peaks,2) % for each type of transition
        
        if isempty(peaks(iPeakFreq).f1)
            disp(['no instances of' peakTypes{iPeakFreq}])
            
        elseif size(peaks(iPeakFreq).f1, 1) < minPkInstances
            disp(['not enough instances of' peakTypes{iPeakFreq}])
            
        else
            meanF1trace = nanmean(peaks(iPeakFreq).f1,1);
            errorF1trace = ste(peaks(iPeakFreq).f1);
            
            meanF2trace = nanmean(peaks(iPeakFreq).f2,1);
            errorF2trace = ste(peaks(iPeakFreq).f2);
            
            figure
            title([parName ' ' runName ' ' peaks(iPeakFreq).type ' peaks'])
            hold on
            mseb([(-peakHalf):(1/512):peakHalf],[meanF1trace; meanF2trace],[errorF1trace; errorF2trace],[],1);
            vline(0,'k','peak')
            xlabel('Time from peak (s)');
            ylabel('Amplitude');
            legend('Low Frequency', 'High Frequency');
            
        end
        
    end
end

