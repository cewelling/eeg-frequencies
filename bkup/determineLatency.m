function [simLatency, latencyError] = determineLatency(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName)

defaultAnalParams %Load default anal params
loadRLSdata %Load RLS data
spectrum = tfr_rls.powspctrm;


for trialIndex = 1:numTrials  %1:numTrials
    
    time = tfr_rls.time;
    
    trialSpect = squeeze(spectrum(trialIndex,:,:,:));
    meanTrialSpect = squeeze(mean(trialSpect,1));

    freq1channels = meanTrialSpect(channelIndices1,:);
    freq2channels = meanTrialSpect(channelIndices2,:);

    for i = 1:size(freq1channels,1) 
        freq1channels(i,:)  = freq1channels(i,:)/nanmax(freq1channels(i,:));
    end
    for i = 1:size(freq2channels,1) 
        freq2channels(i,:)  = freq2channels(i,:)/nanmax(freq2channels(i,:));
    end
    
    realFreqIndices = ~isnan(freq1channels);
    time = time(realFreqIndices == 1);

    %average
    meanfreq1chan = squeeze(mean(freq1channels,1));
    meanfreq2chan = squeeze(mean(freq2channels,1));

    %remove nans
    meanfreq1chan = meanfreq1chan(~isnan(meanfreq1chan));
    meanfreq2chan = meanfreq2chan(~isnan(meanfreq2chan));
    
    for i = smoothWindow+1:(length(meanfreq1chan)-smoothWindow)
        smoothOver = [];
        for j = -smoothWindow:smoothWindow
            smoothOver = [smoothOver meanfreq1chan(i+j)];
        end
        meanfreq1ChanSmoothed(i) = nanmean(smoothOver);
    end
    meanfreq1ChanSmoothed(meanfreq1ChanSmoothed == 0) = NaN;

    for i = smoothWindow+1:(length(meanfreq2chan)-smoothWindow)
        smoothOver = [];
        for j = -smoothWindow:smoothWindow
            smoothOver = [smoothOver meanfreq2chan(i+j)];
        end
        meanfreq2ChanSmoothed(i) = nanmean(smoothOver);
    end
    meanfreq2ChanSmoothed(meanfreq2ChanSmoothed == 0) = NaN;
    
    time = time(1:end-1);

    x = meanfreq1ChanSmoothed(~isnan(meanfreq1ChanSmoothed)); %data to consider
    y = meanfreq2ChanSmoothed(~isnan(meanfreq2ChanSmoothed)); %data to consider
    timepoints = [1:length(x)]; %filter bad timepoints here
    
    x = x(timepoints);
    y = y(timepoints);
    
    time = time(2:end);
    
    %% Compute Cross Correlation - Simulation Run
    
%     Schedule = importdata('../../runScripts/Results/28-Mar-2016_simulationData_stratus01.txt'); % manual
    Schedule = importdata('../../runScripts/Results/25-Apr-2016_simulationData1_stratus11.txt'); % manual
    
    trialSched = Schedule(Schedule(:,1) == (trialIndex),:); %trialIndex+1, 'plus one' because file missing first run of data [for stratus01]
    
    % Determine what trialSched column 6 means
    if trialSched(1,5) == 1 % fast frequency in left eye, so 1 in column 6 means presenting fast frequency
        simulationShow1 = 'fast';
        xMatch = [];
        for i = 1:size(trialSched,1)
            if trialSched(i,6) == 2
                xMatch(i) = 1;
            else
                xMatch(i) = -1;
            end
        end
        yMatch = [];
        for i = 1:size(trialSched,1)
            if trialSched(i,6) == 1
                yMatch(i) = 1;
            else
                yMatch(i) = -1;
            end
        end
    else % fast frequency in right eye, so 1 in column 6 means presenting slow frequency
        simulationShow1 = 'slow';
        xMatch = [];
        for i = 1:size(trialSched,1)
            if trialSched(i,6) == 1
                xMatch(i) = 1;
            else
                xMatch(i) = -1;
            end
        end
        yMatch = [];
        for i = 1:size(trialSched,1)
            if trialSched(i,6) == 2
                yMatch(i) = 1;
            else
                yMatch(i) = -1;
            end
        end
    end
        
    % Turn frequency band power flux into square wave
    
    sqWaveX = x; % slower frequency
    for i = 1:length(sqWaveX)
        if sqWaveX(i) > (noiseThresh*max(sqWaveX))
            sqWaveX(i) = 1;
        else
            sqWaveX(i) = -1;
        end
    end
    
    sqWaveY = y; % faster frequency
    for i = 1:length(sqWaveY)
        if sqWaveY(i) > (noiseThresh*max(sqWaveY))
            sqWaveY(i) = 1;
        else
            sqWaveY(i) = -1;
        end
    end
            
    % Down sample the presentation data
    
    xMatchShort = round(resample(xMatch,length(sqWaveX),length(xMatch)));
    yMatchShort = round(resample(yMatch,length(sqWaveY),length(yMatch)));
    
    % Cross Correlate
    
    % Slow Frequency
        [xcX,lagsX]=xcorr(sqWaveX,xMatchShort);
        [m,i]=max(xcX);
        tauX=lagsX(i);

        % plot signals
        figure
        subplot( 2,1,1)
        plot(time,sqWaveX);
        axis([0 30 -1.25 1.25])
        hold on
        plot(time,xMatchShort, 'r');
        axis([0 30 -1.25 1.25])
        title('Slower Frequency Band Power and Presentation Schedule')
        legend('Freq Band Power','Presentation');
        hold off

        % plot cross-correlation
        subplot( 2, 1,2)
        plot( lagsX, xcX(1:end));
        title('Cross Correlation');
    
    % Fast Frequency
        [xcY,lagsY]=xcorr(sqWaveY,yMatchShort);
        [m,i]=max(xcY);
        tauY=lagsY(i);

        % plot signals
        figure
        subplot( 2,1,1)
        plot(time,sqWaveY);
        axis([0 30 -1.25 1.25])
        hold on
        plot(time,yMatchShort, 'r');
        axis([0 30 -1.25 1.25])
        title('Faster Frequency Band Power and Presentation Schedule')
        legend('Freq Band Power','Presentation');
        hold off

        % plot cross-correlation
        subplot( 2, 1,2)
        plot( lagsY, xcY(1:end));
        title('Cross Correlation');
        
    % Best Latency
        tau = mean([tauX tauY]);
        latency = tau*mean(diff(time)) %number of seconds that frequency band power lags behind presentatin
        
        latencies(trialIndex,:) = latency;
    
    
end

simLatency = mean(latencies);
latencyError = std(latencies);

end