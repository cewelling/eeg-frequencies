% clear all;

defaultAnalParams %Load default anal params
loadRLSdata %Load RLS data

spectrum = tfr_rls.powspctrm;

%% format freq data
omniFreqs = [];
omniDiffs = [];
for trialOfInt = 1:numTrials

    spectrum = tfr_rls.powspctrm;
    trialSpect = squeeze(spectrum(trialOfInt,:,:,:));
    meanTrialSpect = squeeze(mean(trialSpect,1));

    timeInd = tfr_rls.time;
    freqInd = tfr_rls.freq;

    %broadband
    freq1channels = meanTrialSpect(channelIndices1,:); %Freq 1
    freq2channels = meanTrialSpect(channelIndices2,:); %Freq 2

    %normalize: should this be done to the max, median, or mean?
    for i = 1:size(freq1channels,1) 
        freq1channels(i,:)  = freq1channels(i,:)/nanmax(freq1channels(i,:));
    end
    for i = 1:size(freq2channels,1) 
        freq2channels(i,:)  = freq2channels(i,:)/nanmax(freq2channels(i,:));
    end

    %narrowband
    % freq1channels = meanTrialSpect(68,:);
    % freq2channels = meanTrialSpect(85,:);

    buttonPresses = omniData;

    trialPresses = [];
    for i = 1:size(buttonPresses,1)
        if buttonPresses(i,1) == trialOfInt
            trialPresses = [trialPresses; buttonPresses(i,:)];
        end
    end

    %convert epoch start-end into time indices
    lowCount = 1;
    for i = 1:size(trialPresses,1)
        if trialPresses(i,7) == -1
            thisEpoch = find(trialPresses(i,2) < timeInd & timeInd < trialPresses(i,3));
            thisEpoch = find((trialPresses(i,2)-pressLatency) < timeInd & timeInd < (trialPresses(i,3)-pressLatency));
            lowIndices(lowCount,1:2) = [thisEpoch(1) thisEpoch(end)];
            lowCount = lowCount + 1;
        end
    end

    % lowIndices = find(23 < timeInd & timeInd < 24); %4-5
    highCount = 1;
    for i = 1:size(trialPresses,1)
        if trialPresses(i,7) == 1
            thisEpoch = find(trialPresses(i,2) < timeInd & timeInd < trialPresses(i,3));
            thisEpoch = find((trialPresses(i,2)-pressLatency) < timeInd & timeInd < (trialPresses(i,3)-pressLatency));
            highIndices(highCount,1:2) = [thisEpoch(1) thisEpoch(end)];
            highCount = highCount + 1;
        end
    end

    % highIndices = find(20 < timeInd & timeInd < 20.5); %14-15
    %freq channel vals during low/high

    freq1low = nanmean(nanmean(freq1channels(:,lowIndices)));
    freq1high = nanmean(nanmean(freq1channels(:,highIndices)));

    freq2low = nanmean(nanmean(freq2channels(:,lowIndices)));
    freq2high = nanmean(nanmean(freq2channels(:,highIndices)));

%     % plot it
%     figure;
%     bar([freq1low freq1high freq2low freq2high])
%     title(['Averaged across events, trial' num2str(trialOfInt)])


    %% 
    for i = 1:length(lowIndices)
        meanFreq1ChannelsLow(i,:) = nanmean(freq1channels(:,lowIndices(i,1):lowIndices(i,2)),2);
        meanFreq2ChannelsLow(i,:) = nanmean(freq2channels(:,lowIndices(i,1):lowIndices(i,2)),2);
        meanDiffLow(i,:) =  meanFreq1ChannelsLow(i,:) - meanFreq2ChannelsLow(i,:);
    end

    for i = 1:length(highIndices)
        meanFreq1ChannelsHigh(i,:) = nanmean(freq1channels(:,highIndices(i,1):highIndices(i,2)),2);
        meanFreq2ChannelsHigh(i,:) = nanmean(freq2channels(:,highIndices(i,1):highIndices(i,2)),2);
        meanDiffHigh(i,:) =  meanFreq2ChannelsHigh(i,:) - meanFreq1ChannelsHigh(i,:);
    end


    freq1low = nanmean(nanmean(meanFreq1ChannelsLow));
    freq1high = mean(mean(meanFreq1ChannelsHigh));

    freq2low = nanmean(nanmean(meanFreq2ChannelsLow));
    freq2high = nanmean(nanmean(meanFreq2ChannelsHigh));
    
    diff1 = nanmean(nanmean(meanDiffHigh));
    diff2 = nanmean(nanmean(meanDiffLow));

    % plot it
    figure;
    bar([freq1low freq1high freq2low freq2high])
    title(['Averaged across events, trial' num2str(trialOfInt)])

    % plot it
    figure;
    bar([diff1 diff2])
    title(['Difference Score: trial' num2str(trialOfInt)])

    omniFreqs(trialOfInt,:) = [freq1low freq1high freq2low freq2high];
    omniDiffs(trialOfInt,:) = [diff1 diff2];

end

%Plot Average Data: Bars
runSum  = mean(omniFreqs);
stds = std(omniFreqs);
figure
hold on
bar(1:4,runSum)
errorbar(1:4,runSum,stds,'.')
title(['Across Trials'])

%Plot Average Data: Difference Scores
runSum  = mean(omniDiffs);
stds = std(omniDiffs);
figure
hold on
bar(1:2,runSum)
errorbar(1:2,runSum,stds,'.')
title(['Across Trials'])
