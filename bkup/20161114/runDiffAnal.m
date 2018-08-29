function plotFrequencyDifference(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName)

defaultAnalParams %Load default anal params
loadRLSdata %Load RLS data
spectrum = tfr_rls.powspctrm;

%widen the frequency bands
% channelIndices1 = [(channelIndices1-1) channelIndices1 (channelIndices1+1)];
% channelIndices2 = [(channelIndices2-1) channelIndices2 (channelIndices2+1)];

% trimStart = 200;
% spectrum = spectrum(:,:,:,trimStart:end);

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

%     figure;
%     plot(differenceScore1) %plot against tfr_rls.time
%     hold on
%     plot(differenceScore2)



    %% Check against magnitude of the intermodulation frequency
    
    interFreqChannels = meanTrialSpect(interFreqIndices,:);
    
    for i = size(interFreqChannels,1) 
        interFreqChannels(i,:)  = interFreqChannels(i,:)/nanmax(interFreqChannels(i,:));
    end

    meanInterFreqChan = squeeze(mean(interFreqChannels,1));
    meanInterFreqChan = meanInterFreqChan(~isnan(meanInterFreqChan));

    for i = smoothWindow+1:(size(meanInterFreqChan,2) - smoothWindow)
        smoothOver = [];
        for j = -smoothWindow:smoothWindow
            smoothOver = [smoothOver meanInterFreqChan(i+j)];
        end
        meanInterFreqChanSmoothed(i) = nanmean(smoothOver);
    end
    meanInterFreqChanSmoothed(meanInterFreqChanSmoothed == 0) = NaN;

    figure;
    subplot(3,1,1:2);
    plot(time,differenceScore1)
    hold on
    plot(time,differenceScore2)
    box off
    title('Difference Scores')
    subplot(3,1,3);
    plot(time,meanInterFreqChanSmoothed)
    box off
    title('Intermodulation Frequencies')


    %% Bin correlation over time windows and compute correlation

    timeWin = 25; %time bin of sliding window

    x = meanfreq1ChanSmoothed(~isnan(meanfreq1ChanSmoothed)); %data to consider
    y = meanfreq2ChanSmoothed(~isnan(meanfreq2ChanSmoothed)); %data to consider
    timepoints = [1:length(x)]; %filter bad timepoints here
    
    x = x(timepoints);
    y = y(timepoints);
    
    signalX = find(x/max(x) > noiseThresh);
    signalY = find(y/max(y) > noiseThresh);
    
    time = time(2:end);

    keepFreq1 = [];
    keepFreq2 = [];
    for curTime = 1:length(x)-timeWin

        keepFreq1(curTime) = mean(x(curTime:curTime+timeWin));
        keepFreq2(curTime) = mean(y(curTime:curTime+timeWin));

    end
    diffFreq1Freq2 = keepFreq1 - keepFreq2;
    diffFreq2Freq1 = keepFreq2 - keepFreq1;


    %Plot two frequencies across time bins
    figure
    plot(time,x); hold on; plot(time,y); hold on
%     px=[0 size(spectrum,4) size(spectrum,4) 0]; % make closed patch
    px=[0 ceil(max(time)) ceil(max(time)) 0]; % make closed patch
    py=[0 0 noiseThresh*max(x) noiseThresh*max(x)];
    patch(px,py,1,'FaceColor','r','EdgeColor','none');
    alpha(0.25)
    box off
    title(['Signal in each Frequency Channel, Trial ' num2str(trialIndex)])
    legend('Frequency 1','Frequency 2','Noise Threshold');  
    legend('boxoff')
    
    
    
    %WORKING: Alternative summary statistic
    newX = x;
    newY = y;
    [outliersQ3, outliersQ1, modoutliersQ3, modoutliersQ1, IQR,Q] = quartile(newX);  
    sortedX = zeros(1,length(x));
    sortedX(find(x<Q(1))) = -1; sortedX(find(x>Q(3))) = 1;
    
    [outliersQ3, outliersQ1, modoutliersQ3, modoutliersQ1, IQR,Q] = quartile(newY);  
    sortedY = zeros(1,length(y));
    sortedY(find(x<Q(1))) = -1; sortedY(find(x>Q(3))) = 1;
    
    unionZEROSXY = intersect(find(abs(sortedX) == 1), find(abs(sortedY) == 1));
    newX = newX(unionZEROSXY);
    newY = newY(unionZEROSXY);
    
    aboveX = (newX > nanmedian(newX));
    aboveY = (newY > nanmedian(newY));
    sumXY = aboveX+aboveY;
   
    sumXY1(trialIndex) = length(find(sumXY == 1))/length(sumXY);
    
    figure
    plot(aboveX); hold on; plot(aboveY)
    hold on; plot(sumXY)
    
    %When one frequency is high, the other should be zero and visa versa
    figure
    scatter(x,y)

    %Statistics
    [R, t, P] = spear(x',y'); %Are these two timeseries correlated?
    
    diffFreq1Freq2 = x-y;
    diffFreq1Freq2 = diffFreq1Freq2(signalX);
    diffFreq2Freq1 = y-x;
    diffFreq2Freq1 = diffFreq2Freq1(signalY);
    [hDiff,pDiff1] = ttest(diffFreq1Freq2,0) %Is the difference between them typically different from zero
    [hDiff,pDiff2] = ttest(diffFreq2Freq1,0) %Is the difference between them typically different from zero

    %Plot absolute difference from zero across time bins
    meanInterFreqChanSmoothed2 = meanInterFreqChanSmoothed(signalX);
    interModTPs = find(meanInterFreqChanSmoothed2/max(meanInterFreqChanSmoothed2) > threshInterMod);
    
    pad = 2;
    figure
    subplot(1,2,1), bar(1:length(diffFreq1Freq2), diffFreq1Freq2); hold on
    for i = 1:length(interModTPs)
        px=[interModTPs(i)-pad interModTPs(i)+pad interModTPs(i)-pad interModTPs(i)+pad]; % make closed patch
        py=[-1 -1 1 1];
        a = patch(px,py,1,'FaceColor','r','EdgeColor','none');
        alpha(0.2)
        hold on
    end
    box off
    subplot(1,2,2), bar(1:length(diffFreq2Freq1), diffFreq2Freq1); hold on
    for i = 1:length(interModTPs)
        px=[interModTPs(i)-pad interModTPs(i)+pad interModTPs(i)-pad interModTPs(i)+pad]; % make closed patch
        py=[-1 -1 1 1];
        a = patch(px,py,1,'FaceColor','r','EdgeColor','none');
        alpha(0.2)
        hold on
    end
    box off
    legend(a,'Intermodulation Frequencies'); legend('boxoff')
    title(['Difference Scores (Freq 1 - Freq 2) - Noise TPs Excluded, Trial ' num2str(trialIndex)])
    
    omniStats(trialIndex,1:6) = [R P pDiff1 mean(diffFreq1Freq2) pDiff2 mean(diffFreq2Freq1)];

    %% Compute cross Correlation between two frequency time series
%     [c,lags] = xcorr(x,y,'coeff');
%     figure
%     plot(lags,c,'k'), xlabel('lags [steps]'), ylabel('corr')

end

sumXY1

%% Plot Difference Scores Across Trials
data = omniStats(:,[4 6]);
[h,p] = ttest(data,0);
meanData = nanmean(data);
std = ste(data);

figure
bar(1:2,meanData); hold on
errorbar(meanData,std,'b.'); hold on
for i = 1:size(data,2)
    scatter([i*ones(1,length(data(:,i)))]', data(:,i))
end
box off
legend(['p = ' num2str(p)]); legend('boxoff')
title(['Difference Scores Across Trials'])
    
