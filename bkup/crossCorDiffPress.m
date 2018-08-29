function crossCorDiffPress(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName,omniData)

%load the rls data
%calculate difference score across time
%turn difference score into square wave
%load the button press data -> square wave
%cross correlate the reported epochs with the difference score

defaultAnalParams %Load default anal params
loadRLSdata %Load RLS data
spectrum = tfr_rls.powspctrm;

omniStats = NaN*ones(size(spectrum,1),6);
for trialIndex = 1  %1:numTrials
    
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
    time = time(realFreqIndices(1,:) == 1);

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

    differenceScore1 = meanfreq1ChanSmoothed - meanfreq2ChanSmoothed;
    differenceScore2 = meanfreq2ChanSmoothed - meanfreq1ChanSmoothed;
    
    time = time(1:end-1);
    
    %convert to square wave
    for timePt = 1:length(differenceScore1)
        if differenceScore1(timePt) > 0
            differenceScore1(timePt) = 1;
        else
            differenceScore1(timePt) = -1;
        end
        if differenceScore2(timePt) > 0
            differenceScore2(timePt) = 1;
        else
            differenceScore2(timePt) = -1;
        end
    end
    
    %sort and square the button press data
    trialPress = omniData(omniData(:,1) == trialIndex,:);
    
    for i = 1:length(lowFreqReportData.time{1,trial})
        [distance closeEpochRow] = min(abs(meanTimes(meanTimes(:,1) == trial,2)-lowFreqReportData.time{1,trial}(1,i)));
        if lowFreqReportData.time{1,trial}(1,i) > thisButtonReport(closeEpochRow,2) && lowFreqReportData.time{1,1}(1,i) < thisButtonReport(closeEpochRow,3) && thisButtonReport(closeEpochRow,7) == -1
            lowFreqReportData.trial{1,trial}(:,i) = lowFreqReportData.trial{1,trial}(:,i);
        else
            lowFreqReportData.trial{1,trial}(:,i) = 0;
        end
    end
    
    
end
end