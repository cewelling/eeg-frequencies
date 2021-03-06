function [ peaks ] = aggPeaks( parName, runName, date, EEGfile)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% Set-up

% load parameters
analysisParams

% get list of participant's perceptual epochs for this run
[epochs, ~] = getEpochs(parName, runName, date);

% load RLS time amplitude data
[rls_data rls_time] = runRLS(parName, runName, date, EEGfile);

% Preallocate space
fastIndices = [];
slowIndices = [];
fastTrials = [];
slowTrials = [];

%% Find and store dominant state peaks

% Find a reported dominant percept
for iEpoch = 2:(size(epochs,1) - 1) % in first and last epoch, peaks might be cut off
    if epochs(iEpoch, 5) ~= 0
        if epochs(iEpoch, 3) - epochs(iEpoch, 2) > domMin
            
            % Step from the edges of the dom percept & find its peak
            peakIndex = findPeak(epochs(iEpoch, :), rls_data, rls_time);
            
            % Store peaks and the trials in which they appear
            if ~isnan(peakIndex)
                if epochs(iEpoch, 5) == 1 % fast freq dominant
                    fastIndices = [fastIndices peakIndex];
                    fastTrials = [fastTrials epochs(iEpoch, 1)];
                else % slow freq dominant
                    slowIndices = [slowIndices peakIndex];
                    slowTrials = [slowTrials epochs(iEpoch, 1)];
                end
            end
        end
    end
end

peakIndices = {slowIndices fastIndices};
peakTrials = {slowTrials fastTrials};

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

for iTrial = 1:size(rls_data(1).amp, 2)   
     for iFreq = 1:length(rls_data);
         
         RLSamp = (rls_data(iFreq).amp{iTrial})';
        
%         % take normalizing/standardizing step
%         if strcmp(normType, 'mean')
             nRLSamp = RLSamp - mean(RLSamp);
%         elseif strcmp(normType, 'norm')
%             normRLSamp = 2*(RLSamp - min(RLSamp))/(max(RLSamp) - min(RLSamp)) - 1;
%         elseif strcmp(normType, 'z')
%             normRLSamp = zscore(RLSamp);
%         elseif strcmp(normType, 'none')
%             normRLSamp = RLSamp;
%         end
%         
         nRLS_data(iFreq).amp{iTrial} = nRLSamp';        
     end
%     
 end

% nRLS_data = rls_data;

%% Average peaks

halfLength = 512*peakHalf; % Remember we have 512 data points per second
fullLength = halfLength*2+1; % twice the half length plus the transition time
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
        
        % epoch is cut off by beginning, pad beginning with NaNs
        if (peakIndex - halfLength) < 1
            epochRLSf1 = RLSf1(1:peakIndex + halfLength);
            epochRLSf2 = RLSf2(1:peakIndex + halfLength);
            padSize = fullLength - length(epochRLSf1);

            epochRLSf1 = [nan(1,padSize), epochRLSf1];
            epochRLSf2 = [nan(1,padSize), epochRLSf2];
            
        % epoch is cut off by end, pad end with NaNs   
        elseif (peakIndex + halfLength) > length(RLSf1) - (cutoffWindow / 1000 * 512); % remember end is padded
            epochRLSf1 = RLSf1(peakIndex - halfLength : length(RLSf1) - (cutoffWindow / 1000 * 512));
            epochRLSf2 = RLSf2(peakIndex - halfLength : length(RLSf2) - (cutoffWindow / 1000 * 512));
            padSize = fullLength - length(epochRLSf1);
            
            epochRLSf1 = [epochRLSf1,nan(1,padSize)];          
            epochRLSf2 = [epochRLSf2,nan(1,padSize)];
            
        % epoch is not cut off    
        else            
            epochRLSf1 = RLSf1(peakIndex-halfLength:peakIndex + halfLength);
            epochRLSf2 = RLSf2(peakIndex-halfLength:peakIndex + halfLength);
        end
        
            f1Epochs = [f1Epochs; epochRLSf1];
            f2Epochs = [f2Epochs; epochRLSf2];
            
    end
    
    peaks(iPeakFreq).type = peakTypes{iPeakFreq};
    peaks(iPeakFreq).f1 = f1Epochs;
    peaks(iPeakFreq).f2 = f2Epochs;
    
end

%% Plot averaged peaks

if strcmp(avgPksPlotOrNot, 'yes')
    for iPeakFreq = 1:size(peaks,2) % for each type of transition
        
        if isempty(peaks(iPeakFreq).f1)
            disp(['no instances of' peaks(iFreqPeak.type)])
            
        elseif size(peaks(iPeakFreq).f1, 1) < minPkInstances
            disp(['not enough instances of' peaks(iFreqPeak.type)])
            
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

